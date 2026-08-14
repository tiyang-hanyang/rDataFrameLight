import argparse
import warnings
from pathlib import Path

import awkward as ak
import numpy as np
import torch
import uproot

from feature_grouping import select_active_features
from predictSignalRootToNpy import build_model, extract_features_with_ids, run_inference
from runSignalRootInferenceFromSampleJson import (
    build_jobs,
    sanitize_namespace,
)
from writeSignalPredictionToRoot import DEFAULT_FLOAT, DEFAULT_INT, DEFAULT_STRING


COMMON_JES_SOURCES = [
    "CMS_scale_j_FlavorQCD",
    "CMS_scale_j_RelativeBal",
    "CMS_scale_j_HF",
    "CMS_scale_j_BBEC1",
    "CMS_scale_j_EC2",
    "CMS_scale_j_Absolute",
]

YEAR_JES_SOURCES = {
    "2022": [
        "CMS_scale_j_Absolute_2022",
        "CMS_scale_j_HF_2022",
        "CMS_scale_j_EC2_2022",
        "CMS_scale_j_RelativeSample_2022",
        "CMS_scale_j_BBEC1_2022",
    ],
    "2022EE": [
        "CMS_scale_j_Absolute_2022EE",
        "CMS_scale_j_HF_2022EE",
        "CMS_scale_j_EC2_2022EE",
        "CMS_scale_j_RelativeSample_2022EE",
        "CMS_scale_j_BBEC1_2022EE",
    ],
    "2023": [
        "CMS_scale_j_Absolute_2023",
        "CMS_scale_j_HF_2023",
        "CMS_scale_j_EC2_2023",
        "CMS_scale_j_RelativeSample_2023",
        "CMS_scale_j_BBEC1_2023",
    ],
    "2023BPix": [
        "CMS_scale_j_Absolute_2023BPix",
        "CMS_scale_j_HF_2023BPix",
        "CMS_scale_j_EC2_2023BPix",
        "CMS_scale_j_RelativeSample_2023BPix",
        "CMS_scale_j_BBEC1_2023BPix",
    ],
    "2024": [
        "CMS_scale_j_Absolute_2024",
        "CMS_scale_j_HF_2024",
        "CMS_scale_j_EC2_2024",
        "CMS_scale_j_RelativeSample_2024",
        "CMS_scale_j_BBEC1_2024",
    ],
}

DEFAULT_BTAG_THRESHOLDS = {
    "2022": 0.245,
    "2022EE": 0.2605,
    "2023": 0.1917,
    "2023BPix": 0.1919,
    "2024": 0.1272,
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run signal DNN inference once per ROOT file and write nominal plus JME systematic "
            "score branches to one NanoAOD-like output ROOT. Data receives nominal only."
        )
    )
    parser.add_argument("--sample-json", "--sampleJson", dest="sample_json", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--tree-name", default="Events")
    parser.add_argument("--prefix", default="Signal")
    parser.add_argument("--campaign", choices=sorted(YEAR_JES_SOURCES), required=True)
    parser.add_argument("--sample-keys", nargs="+", default=None)
    parser.add_argument("--batch-size", type=int, default=2048)
    parser.add_argument("--device", default=None)
    parser.add_argument("--step-size", default="100 MB")
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--stop-on-error", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--btag-threshold", type=float, default=None)
    parser.add_argument("--jet-pt-threshold", type=float, default=30.0)
    parser.add_argument("--scale-jet-mass-with-pt", action="store_true")
    parser.add_argument(
        "--include-jes",
        action="store_true",
        default=True,
        help="Include JES score variations for MC. Enabled by default.",
    )
    parser.add_argument(
        "--no-jes",
        dest="include_jes",
        action="store_false",
        help="Disable JES score variations.",
    )
    parser.add_argument(
        "--include-jer",
        action="store_true",
        default=True,
        help="Include JER score variations for MC. Enabled by default.",
    )
    parser.add_argument(
        "--no-jer",
        dest="include_jer",
        action="store_false",
        help="Disable JER score variations.",
    )
    parser.add_argument(
        "--include-muon-scale",
        action="store_true",
        default=True,
        help="Include muon Rochester scale score variations for MC. Enabled by default.",
    )
    parser.add_argument(
        "--no-muon-scale",
        dest="include_muon_scale",
        action="store_false",
        help="Disable muon Rochester scale score variations.",
    )
    parser.add_argument(
        "--include-muon-resol",
        action="store_true",
        default=True,
        help="Include muon Rochester resolution score variations for MC. Enabled by default.",
    )
    parser.add_argument(
        "--no-muon-resol",
        dest="include_muon_resol",
        action="store_false",
        help="Disable muon Rochester resolution score variations.",
    )
    return parser.parse_args()


def resolve_device(device_name):
    if device_name is not None:
        return torch.device(device_name)
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


def is_data_job(job):
    return any(part == "data" for part in Path(job["relative_dir"]).parts)


def build_variations(args, is_data):
    variations = [("nominal", args.prefix, {})]
    if is_data:
        return variations

    if args.include_jes:
        for source in COMMON_JES_SOURCES + YEAR_JES_SOURCES[args.campaign]:
            for direction in ("up", "down"):
                label = f"{source}_{direction}"
                variations.append(
                    (
                        label,
                        f"{args.prefix}_{label}",
                        {
                            "jet_pt_shift_branch": source,
                            "jet_pt_shift_direction": direction,
                        },
                    )
                )

    if args.include_jer:
        for direction in ("up", "down"):
            label = f"JER_corr_{direction}"
            variations.append(
                (
                    label,
                    f"{args.prefix}_{label}",
                    {
                        "jer_direction": direction,
                    },
                )
            )

    if args.include_muon_scale:
        for direction, branch_name in (
            ("up", "Muon_pt_Rscale_up"),
            ("down", "Muon_pt_Rscale_dn"),
        ):
            label = f"MuonRoch_scale_{direction}"
            variations.append(
                (
                    label,
                    f"{args.prefix}_{label}",
                    {
                        "muon_pt_branch": branch_name,
                    },
                )
            )

    if args.include_muon_resol:
        for direction, branch_name in (
            ("up", "Muon_pt_Rcorr_resolup"),
            ("down", "Muon_pt_Rcorr_resoldn"),
        ):
            label = f"MuonRoch_resol_{direction}"
            variations.append(
                (
                    label,
                    f"{args.prefix}_{label}",
                    {
                        "muon_pt_branch": branch_name,
                    },
                )
            )
    return variations


def run_prediction_for_variation(args, model, contract, device, input_path, tree_name, variation_kwargs):
    ids, features = extract_features_with_ids(
        input_path,
        tree_name,
        jet_pt_shift_branch=variation_kwargs.get("jet_pt_shift_branch"),
        jet_pt_shift_direction=variation_kwargs.get("jet_pt_shift_direction", "up"),
        jer_direction=variation_kwargs.get("jer_direction"),
        muon_pt_branch=variation_kwargs.get("muon_pt_branch"),
        jet_pt_threshold=args.jet_pt_threshold,
        btag_threshold=args.btag_threshold,
        scale_jet_mass_with_pt=args.scale_jet_mass_with_pt,
    )
    if contract["feature_layout"] == "compact":
        features = select_active_features(features, contract["active_feature_groups"])
    class_names = list(contract["task_label_names"])
    scores = run_inference(model, features, args.batch_size, device, len(class_names))
    return ids, scores, class_names


def build_prediction_table(ids, scores, class_names):
    table = {}
    for index, key in enumerate(ids):
        table[(int(key[0]), int(key[1]), int(key[2]))] = scores[index]
    return table


def build_score_columns(arrays, prediction_tables, class_names_by_prefix):
    runs = ak.to_numpy(arrays["run"])
    lumis = ak.to_numpy(arrays["luminosityBlock"])
    events = ak.to_numpy(arrays["event"])
    n_events = len(runs)
    columns = {}
    matched_counts = {}

    for prefix, table in prediction_tables.items():
        class_names = class_names_by_prefix[prefix]
        prediction = np.full(n_events, DEFAULT_INT, dtype=np.int32)
        best_score = np.full(n_events, DEFAULT_FLOAT, dtype=np.float32)
        max_name_length = max([len(DEFAULT_STRING)] + [len(name) for name in class_names])
        predicted_class_name = np.full(n_events, DEFAULT_STRING, dtype=f"<U{max_name_length}")
        scores = {class_name: np.full(n_events, DEFAULT_FLOAT, dtype=np.float32) for class_name in class_names}
        matched = 0

        for index, key in enumerate(zip(runs, lumis, events)):
            score_row = table.get((int(key[0]), int(key[1]), int(key[2])))
            if score_row is None:
                continue
            matched += 1
            pred = int(np.argmax(score_row))
            prediction[index] = pred
            best_score[index] = np.float32(score_row[pred])
            predicted_class_name[index] = class_names[pred]
            for class_index, class_name in enumerate(class_names):
                scores[class_name][index] = np.float32(score_row[class_index])

        columns[f"{prefix}_prediction"] = prediction
        columns[f"{prefix}_bestScore"] = best_score
        columns[f"{prefix}_predictedClassName"] = ak.Array(predicted_class_name.tolist())
        columns[f"{prefix}_isMatched"] = prediction >= 0
        for class_name in class_names:
            columns[f"{prefix}_score_{class_name}"] = scores[class_name]
        matched_counts[prefix] = matched
    return columns, matched_counts


def write_scored_root(input_path, output_path, tree_name, prediction_tables, class_names_by_prefix, step_size):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    total_input_events = 0
    total_matched = {prefix: 0 for prefix in prediction_tables}
    has_written = False

    with uproot.open(input_path) as input_file, uproot.recreate(output_path) as output_file:
        for key in input_file.keys(cycle=False):
            if key == tree_name:
                continue
            output_file[key] = input_file[key]

        for arrays in uproot.iterate(f"{input_path}:{tree_name}", library="ak", step_size=step_size):
            total_input_events += len(arrays["run"])
            new_columns, matched_counts = build_score_columns(arrays, prediction_tables, class_names_by_prefix)
            for prefix, matched in matched_counts.items():
                total_matched[prefix] += matched

            merged = {field: arrays[field] for field in arrays.fields}
            merged.update(new_columns)
            if not has_written:
                output_file[tree_name] = merged
                has_written = True
            else:
                output_file[tree_name].extend(merged)

    return total_input_events, total_matched


def main():
    args = parse_args()
    if args.btag_threshold is None:
        args.btag_threshold = DEFAULT_BTAG_THRESHOLDS[args.campaign]
    device = resolve_device(args.device)
    model, contract = build_model(args.model, device)

    output_root = Path(args.output_dir)
    jobs = build_jobs(args.sample_json, args.sample_keys, args.tree_name)
    if not jobs:
        warnings.warn(f"No valid files were found in {args.sample_json}. Nothing to do.", stacklevel=2)
        return

    failed_jobs = []
    for job in jobs:
        relative_dir = Path(job["relative_dir"])
        input_path = job["input_path"]
        output_dir = output_root / relative_dir
        output_path = output_dir / job["file_name"]
        data_job = is_data_job(job)
        variations = build_variations(args, data_job)

        print(f"sample_key: {job['sample_key']}")
        print(f"input: {input_path}")
        print(f"output: {output_path}")
        print(f"data: {data_job}")
        print(f"variations: {len(variations)}")

        if args.skip_existing and output_path.is_file():
            print("skip: output already exists")
            continue
        if args.dry_run:
            continue

        try:
            prediction_tables = {}
            class_names_by_prefix = {}
            for label, prefix, kwargs in variations:
                print(f"  infer: {label} -> {prefix}")
                ids, scores, class_names = run_prediction_for_variation(
                    args,
                    model,
                    contract,
                    device,
                    input_path,
                    args.tree_name,
                    kwargs,
                )
                prediction_tables[prefix] = build_prediction_table(ids, scores, class_names)
                class_names_by_prefix[prefix] = class_names

            total_events, total_matched = write_scored_root(
                input_path,
                output_path,
                args.tree_name,
                prediction_tables,
                class_names_by_prefix,
                args.step_size,
            )
            print(f"  input events written: {total_events}")
            for prefix, matched in total_matched.items():
                print(f"  matched[{prefix}]: {matched}")
        except Exception as exc:
            if args.stop_on_error:
                raise
            warnings.warn(f"Failed for {input_path}: {exc}", stacklevel=2)
            failed_jobs.append((input_path, str(exc)))

    print(f"output root: {output_root}")
    print(f"failed files: {len(failed_jobs)}")
    for input_path, reason in failed_jobs:
        print(f"failed: {input_path} :: {reason}")


if __name__ == "__main__":
    main()
