import argparse
import json
import subprocess
import sys
import warnings
from pathlib import Path, PurePosixPath

import uproot


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Batch wrapper around predictSignalRootToNpy.py and writeSignalPredictionToRoot.py "
            "using a sample json with dir/file maps."
        )
    )
    parser.add_argument("--sample-json", "--sampleJson", dest="sample_json", required=True)
    parser.add_argument("--model", required=True, help="Trained signal model checkpoint.")
    parser.add_argument("--output-dir", required=True, help="Root output directory for predictions and scored ROOT files.")
    parser.add_argument("--tree-name", default="Events")
    parser.add_argument("--prefix", default="Signal", help="Prefix for branches written to ROOT.")
    parser.add_argument(
        "--namespace",
        default=None,
        help=(
            "Optional namespace used for internal prediction-cache files. Defaults to the sample json stem. "
            "Scored ROOT files are written directly under output-dir/data/... or output-dir/mc/...."
        ),
    )
    parser.add_argument(
        "--sample-keys",
        nargs="+",
        default=None,
        help="Optional subset of sample keys from the sample json to run.",
    )
    parser.add_argument(
        "--python-bin",
        default=sys.executable if sys.executable else "python3",
        help="Python executable used to invoke the underlying scripts.",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip predict/write steps when their expected outputs already exist.",
    )
    parser.add_argument(
        "--skip-predict",
        action="store_true",
        help="Skip ROOT -> prediction.npy and require existing prediction outputs.",
    )
    parser.add_argument(
        "--skip-write",
        action="store_true",
        help="Skip writing scored ROOT files.",
    )
    parser.add_argument(
        "--step-size",
        default="100 MB",
        help="Chunk size forwarded to writeSignalPredictionToRoot.py.",
    )
    parser.add_argument("--batch-size", type=int, default=2048, help="Inference batch size.")
    parser.add_argument("--device", default=None, help="cpu or cuda")
    parser.add_argument("--jet-pt-shift-branch", default=None, help="Forwarded to predictSignalRootToNpy.py.")
    parser.add_argument("--jet-pt-shift-direction", choices=["up", "down"], default="up")
    parser.add_argument("--jer-direction", choices=["up", "down"], default=None)
    parser.add_argument("--muon-pt-branch", default=None, help="Forwarded to predictSignalRootToNpy.py.")
    parser.add_argument("--jet-pt-threshold", type=float, default=30.0)
    parser.add_argument("--btag-threshold", type=float, default=None)
    parser.add_argument("--scale-jet-mass-with-pt", action="store_true")
    parser.add_argument(
        "--stop-on-error",
        action="store_true",
        help="Stop immediately on the first failed input file. Default behavior is to warn and continue.",
    )
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def load_sample_json(sample_json_path):
    with open(sample_json_path, "r", encoding="utf-8") as handle:
        payload = json.load(handle)
    dir_map = payload.get("dir")
    file_map = payload.get("file")
    if not isinstance(dir_map, dict) or not isinstance(file_map, dict):
        raise RuntimeError(
            f"Sample json must contain object fields 'dir' and 'file': {sample_json_path}"
        )
    return dir_map, file_map


def join_sample_path(base_dir, file_name):
    base_dir = str(base_dir)
    if "/" in base_dir:
        return str(PurePosixPath(base_dir) / file_name)
    return str(Path(base_dir) / file_name)


def infer_relative_sample_dir(sample_dir, sample_key):
    normalized = str(sample_dir).replace("\\", "/").rstrip("/")
    parts = [part for part in normalized.split("/") if part]
    for marker in ("data", "mc"):
        if marker in parts:
            marker_index = parts.index(marker)
            return Path(*parts[marker_index:])
    return Path(sample_key)


def sanitize_namespace(name):
    cleaned = "".join(char if char.isalnum() or char in ("-", "_", ".") else "_" for char in str(name))
    return cleaned.strip("._") or "inference_batch"


def inspect_root_file(input_path, requested_tree_name):
    with uproot.open(input_path) as root_file:
        keys = set(root_file.keys(cycle=False))
        has_gen_weight_sum = "genWeightSum" in keys
        has_requested_tree = requested_tree_name in keys
        return has_requested_tree, has_gen_weight_sum


def build_jobs(sample_json_path, selected_sample_keys, requested_tree_name):
    dir_map, file_map = load_sample_json(sample_json_path)
    available_keys = sorted(set(dir_map) | set(file_map))
    if selected_sample_keys is None:
        chosen_keys = available_keys
    else:
        chosen_keys = []
        for sample_key in selected_sample_keys:
            if sample_key not in available_keys:
                raise RuntimeError(
                    f"Sample key '{sample_key}' not found in {sample_json_path}. "
                    f"Available keys: {', '.join(available_keys)}"
                )
            chosen_keys.append(sample_key)

    jobs = []
    for sample_key in chosen_keys:
        sample_dir = dir_map.get(sample_key)
        file_list = file_map.get(sample_key, [])
        if not sample_dir:
            warnings.warn(
                f"Skipping sample key '{sample_key}' because 'dir' is missing in {sample_json_path}",
                stacklevel=2,
            )
            continue
        if not file_list:
            warnings.warn(
                f"Skipping sample key '{sample_key}' because 'file' is missing or empty in {sample_json_path}",
                stacklevel=2,
            )
            continue
        relative_dir = infer_relative_sample_dir(sample_dir, sample_key)
        for file_name in file_list:
            input_path = join_sample_path(sample_dir, file_name)
            try:
                has_requested_tree, has_gen_weight_sum = inspect_root_file(input_path, requested_tree_name)
            except Exception as exc:
                warnings.warn(
                    f"Skipping file '{input_path}' because it could not be opened: {exc}",
                    stacklevel=2,
                )
                continue
            if not has_requested_tree:
                warnings.warn(
                    f"Skipping file '{input_path}' because it is missing required object(s): {requested_tree_name}",
                    stacklevel=2,
                )
                continue
            if not has_gen_weight_sum:
                warnings.warn(
                    f"File '{input_path}' does not contain genWeightSum. Continuing without it.",
                    stacklevel=2,
                )
            jobs.append(
                {
                    "sample_key": sample_key,
                    "input_path": input_path,
                    "relative_dir": relative_dir,
                    "file_name": file_name,
                }
            )
    return jobs


def build_predict_command(args, predictor_path, input_path, prediction_path):
    command = [
        args.python_bin,
        str(predictor_path),
        "--input",
        input_path,
        "--model",
        args.model,
        "--output",
        str(prediction_path),
        "--tree-name",
        args.tree_name,
        "--batch-size",
        str(args.batch_size),
    ]
    if args.device is not None:
        command.extend(["--device", args.device])
    if args.jet_pt_shift_branch is not None:
        command.extend(["--jet-pt-shift-branch", args.jet_pt_shift_branch])
        command.extend(["--jet-pt-shift-direction", args.jet_pt_shift_direction])
    if args.jer_direction is not None:
        command.extend(["--jer-direction", args.jer_direction])
    if args.muon_pt_branch is not None:
        command.extend(["--muon-pt-branch", args.muon_pt_branch])
    command.extend(["--jet-pt-threshold", str(args.jet_pt_threshold)])
    if args.btag_threshold is not None:
        command.extend(["--btag-threshold", str(args.btag_threshold)])
    if args.scale_jet_mass_with_pt:
        command.append("--scale-jet-mass-with-pt")
    return command


def build_write_command(args, writer_path, input_path, prediction_path, output_root_path):
    command = [
        args.python_bin,
        str(writer_path),
        "--input",
        input_path,
        "--prediction",
        str(prediction_path),
        "--output",
        str(output_root_path),
        "--tree-name",
        args.tree_name,
        "--prefix",
        args.prefix,
        "--step-size",
        args.step_size,
    ]
    return command


def run_command(command, dry_run):
    print("cmd:", " ".join(command))
    if not dry_run:
        subprocess.run(command, check=True)


def require_existing_file(path, context):
    if not Path(path).is_file():
        raise RuntimeError(f"Missing {context}: {path}")


def main():
    args = parse_args()
    script_dir = Path(__file__).resolve().parent
    predictor_path = script_dir / "predictSignalRootToNpy.py"
    writer_path = script_dir / "writeSignalPredictionToRoot.py"

    sample_json_path = Path(args.sample_json)
    namespace = sanitize_namespace(args.namespace or sample_json_path.stem)
    output_root = Path(args.output_dir)
    prediction_root = output_root / "_predictions" / namespace
    scored_root = output_root
    prediction_root.mkdir(parents=True, exist_ok=True)
    scored_root.mkdir(parents=True, exist_ok=True)

    jobs = build_jobs(args.sample_json, args.sample_keys, args.tree_name)
    if not jobs:
        warnings.warn(
            f"No valid files were found in {args.sample_json}. Nothing to do.",
            stacklevel=2,
        )
        return

    failed_jobs = []
    for job in jobs:
        relative_dir = Path(job["relative_dir"])
        input_path = job["input_path"]
        file_name = job["file_name"]
        file_stem = Path(file_name).stem

        prediction_dir = prediction_root / relative_dir
        output_root_dir = output_root / relative_dir
        prediction_dir.mkdir(parents=True, exist_ok=True)
        output_root_dir.mkdir(parents=True, exist_ok=True)

        prediction_path = prediction_dir / f"{file_stem}.npy"
        output_root_path = output_root_dir / file_name

        print(f"sample_key: {job['sample_key']}")
        print(f"input: {input_path}")
        print(f"prediction: {prediction_path}")
        print(f"scored_root: {output_root_path}")

        if args.skip_predict:
            require_existing_file(prediction_path, f"prediction file for {input_path}")
        else:
            if args.skip_existing and prediction_path.is_file():
                print("skip predict: prediction output already exists")
            else:
                try:
                    run_command(
                        build_predict_command(args, predictor_path, input_path, prediction_path),
                        args.dry_run,
                    )
                except subprocess.CalledProcessError as exc:
                    if args.stop_on_error:
                        raise
                    warnings.warn(
                        f"Predict failed for {input_path}: {exc}. Skipping this file.",
                        stacklevel=2,
                    )
                    failed_jobs.append((input_path, "predict"))
                    continue
                except Exception as exc:
                    if args.stop_on_error:
                        raise
                    warnings.warn(
                        f"Predict failed for {input_path}: {exc}. Skipping this file.",
                        stacklevel=2,
                    )
                    failed_jobs.append((input_path, "predict"))
                    continue

        if args.skip_write:
            continue

        if args.skip_existing and output_root_path.is_file():
            print("skip write: scored ROOT already exists")
            continue

        if not args.dry_run:
            try:
                require_existing_file(prediction_path, f"prediction file for {input_path}")
            except Exception as exc:
                if args.stop_on_error:
                    raise
                warnings.warn(
                    f"Write skipped for {input_path}: {exc}",
                    stacklevel=2,
                )
                failed_jobs.append((input_path, "write_missing_prediction"))
                continue
        try:
            run_command(
                build_write_command(args, writer_path, input_path, prediction_path, output_root_path),
                args.dry_run,
            )
        except subprocess.CalledProcessError as exc:
            if args.stop_on_error:
                raise
            warnings.warn(
                f"Write failed for {input_path}: {exc}. Skipping this file.",
                stacklevel=2,
            )
            failed_jobs.append((input_path, "write"))
            continue
        except Exception as exc:
            if args.stop_on_error:
                raise
            warnings.warn(
                f"Write failed for {input_path}: {exc}. Skipping this file.",
                stacklevel=2,
            )
            failed_jobs.append((input_path, "write"))
            continue

    print(f"namespace: {namespace}")
    print(f"prediction root: {prediction_root}")
    print(f"scored ROOT root: {output_root}")
    print(f"failed files: {len(failed_jobs)}")
    for input_path, stage in failed_jobs:
        print(f"failed[{stage}]: {input_path}")


if __name__ == "__main__":
    main()
