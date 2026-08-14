import argparse
import json
import subprocess
import sys
import warnings
from pathlib import Path, PurePosixPath


DEFAULT_CHANNEL_SPECS = [
    {"channel": "TTHH_DL", "sample_keys": ["TTHH_2B2W_DL", "TTHH_DL_2B2W_batch1"]},
    {"channel": "TTHH_SL", "sample_keys": ["TTHH_2B2W_SL", "TTHH_SL_2B2W_batch1"]},
    {"channel": "ttbarDL", "sample_keys": ["ttbarDL"]},
    {"channel": "ttbarSL", "sample_keys": ["ttbarSL"]},
    {"channel": "TTBB_DL", "sample_keys": ["TTBB_DL"]},
    {"channel": "TTBB_SL", "sample_keys": ["TTBB_SL"]},
    {"channel": "TTHBB", "sample_keys": ["TTHBB"]},
    {"channel": "TTW", "sample_keys": ["TTW"]},
    {"channel": "TTZ_high", "sample_keys": ["TTZ_high"]},
    {"channel": "TTZ_low", "sample_keys": ["TTZ_low"]},
    {"channel": "TTHnonBB", "sample_keys": ["TTHnonBB"]},
    {"channel": "TTTT", "sample_keys": ["TTTT"]},
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Batch wrapper around signalExtractionInputFromNanoAOD.py using a sample json."
    )
    parser.add_argument("--sample-json", "--sampleJson", dest="sample_json", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--era", required=True)
    parser.add_argument("--tree-name", default="Events")
    parser.add_argument("--verbose", type=int, default=1)
    parser.add_argument(
        "--channels",
        nargs="+",
        default=None,
        help="Subset of configured channel names to run.",
    )
    parser.add_argument(
        "--python-bin",
        default=sys.executable if sys.executable else "python3",
        help="Python executable used to invoke signalExtractionInputFromNanoAOD.py",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing them.",
    )
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


def resolve_channel_specs(selected_channels):
    spec_map = {item["channel"]: item for item in DEFAULT_CHANNEL_SPECS}
    if selected_channels is None:
        return DEFAULT_CHANNEL_SPECS
    resolved = []
    for channel_name in selected_channels:
        if channel_name not in spec_map:
            raise RuntimeError(f"Unsupported channel: {channel_name}")
        resolved.append(spec_map[channel_name])
    return resolved


def build_file_jobs(sample_json_path, selected_channels):
    dir_map, file_map = load_sample_json(sample_json_path)
    jobs = []
    for spec in resolve_channel_specs(selected_channels):
        file_jobs = []
        for sample_key in spec["sample_keys"]:
            input_dir = dir_map.get(sample_key)
            input_files = file_map.get(sample_key, [])
            if not input_dir:
                warnings.warn(
                    f"Skipping sample key '{sample_key}' because 'dir' is missing in {sample_json_path}",
                    stacklevel=2,
                )
                continue
            if not input_files:
                warnings.warn(
                    f"Skipping sample key '{sample_key}' because 'file' is missing or empty in {sample_json_path}",
                    stacklevel=2,
                )
                continue
            for input_file in input_files:
                file_jobs.append(
                    {
                        "channel": spec["channel"],
                        "sample_key": sample_key,
                        "input_path": join_sample_path(input_dir, input_file),
                    }
                )
        if not file_jobs:
            warnings.warn(
                f"Skipping channel '{spec['channel']}' because no matching inputs were found in {sample_json_path}",
                stacklevel=2,
            )
            continue
        jobs.append({"channel": spec["channel"], "files": file_jobs})
    return jobs


def build_command(extractor_path, python_bin, input_path, output_path, sample_name, args):
    return [
        python_bin,
        str(extractor_path),
        "--input",
        input_path,
        "--output",
        str(output_path),
        "--sample-name",
        sample_name,
        "--era",
        args.era,
        "--tree-name",
        args.tree_name,
        "--verbose",
        str(args.verbose),
    ]


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    script_dir = Path(__file__).resolve().parent
    extractor_path = script_dir / "signalExtractionInputFromNanoAOD.py"
    jobs = build_file_jobs(args.sample_json, args.channels)
    if not jobs:
        warnings.warn(
            f"No extractable channels were found in {args.sample_json}. Nothing to do.",
            stacklevel=2,
        )
        return

    for job in jobs:
        channel_dir = output_dir / job["channel"]
        channel_dir.mkdir(parents=True, exist_ok=True)
        print(f"channel: {job['channel']}")
        for file_job in job["files"]:
            input_path = file_job["input_path"]
            sample_key = file_job["sample_key"]
            output_path = channel_dir / f"{sample_key}.npy"
            command = build_command(
                extractor_path=extractor_path,
                python_bin=args.python_bin,
                input_path=input_path,
                output_path=output_path,
                sample_name=sample_key,
                args=args,
            )
            print("cmd:", " ".join(command))
            if not args.dry_run:
                subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
