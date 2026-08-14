import argparse
from pathlib import Path

import awkward as ak
import numpy as np
import uproot


DEFAULT_INT = -1
DEFAULT_FLOAT = -999.0
SCORE_BRANCHES = ["b1b2", "b1b3", "b1b4", "b2b3", "b2b4", "b3b4"]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Append ttbb BJA prediction branches from a structured .npy file to a new ROOT file."
    )
    parser.add_argument("--input", required=True, help="Input ROOT file.")
    parser.add_argument("--prediction", required=True, help="Structured .npy produced by predictBJARootToNpy.py.")
    parser.add_argument("--output", required=True, help="Output ROOT file.")
    parser.add_argument("--tree-name", default="Events", help="Input/output tree name.")
    parser.add_argument("--prefix", default="BJA", help="Prefix for appended branches.")
    parser.add_argument("--step-size", default="100 MB", help="Chunk size passed to uproot.iterate.")
    parser.add_argument(
        "--skip-copy-objects",
        action="store_true",
        help="Do not copy non-tree top-level ROOT objects such as genWeightSum.",
    )
    return parser.parse_args()


def load_prediction_table(prediction_path):
    prediction = np.load(prediction_path, allow_pickle=False)
    required_fields = {"run", "luminosityBlock", "event", "prediction", *SCORE_BRANCHES}
    missing = required_fields.difference(prediction.dtype.names or ())
    if missing:
        raise RuntimeError("Prediction file is missing required fields: " + ", ".join(sorted(missing)))

    table = {}
    for row in prediction:
        key = (int(row["run"]), int(row["luminosityBlock"]), int(row["event"]))
        table[key] = row
    return table


def build_new_columns(arrays, table, prefix):
    runs = ak.to_numpy(arrays["run"])
    lumis = ak.to_numpy(arrays["luminosityBlock"])
    events = ak.to_numpy(arrays["event"])
    n_events = len(runs)

    prediction = np.full(n_events, DEFAULT_INT, dtype=np.int32)
    scores = {branch: np.full(n_events, DEFAULT_FLOAT, dtype=np.float32) for branch in SCORE_BRANCHES}

    matched = 0
    for index, key in enumerate(zip(runs, lumis, events)):
        row = table.get((int(key[0]), int(key[1]), int(key[2])))
        if row is None:
            continue
        matched += 1
        prediction[index] = int(row["prediction"])
        for branch in SCORE_BRANCHES:
            scores[branch][index] = np.float32(row[branch])

    new_columns = {
        f"{prefix}_prediction": prediction,
        f"{prefix}_isMatched": prediction >= 0,
    }
    for branch in SCORE_BRANCHES:
        new_columns[f"{prefix}_{branch}"] = scores[branch]
    return new_columns, matched


def main():
    args = parse_args()

    input_path = Path(args.input)
    prediction_path = Path(args.prediction)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    table = load_prediction_table(prediction_path)
    total_input_events = 0
    total_matched_events = 0
    has_written = False

    with uproot.open(input_path) as input_file, uproot.recreate(output_path) as output_file:
        if not args.skip_copy_objects:
            for key in input_file.keys(cycle=False):
                if key == args.tree_name:
                    continue
                output_file[key] = input_file[key]

        for arrays in uproot.iterate(f"{input_path}:{args.tree_name}", library="ak", step_size=args.step_size):
            total_input_events += len(arrays["run"])
            new_columns, matched = build_new_columns(arrays, table, args.prefix)
            total_matched_events += matched

            merged = {field: arrays[field] for field in arrays.fields}
            merged.update(new_columns)

            if not has_written:
                output_file[args.tree_name] = merged
                has_written = True
            else:
                output_file[args.tree_name].extend(merged)

    print(f"prediction rows loaded: {len(table)}")
    print(f"input events written: {total_input_events}")
    print(f"matched events filled: {total_matched_events}")
    print(f"copied extra ROOT objects: {'no' if args.skip_copy_objects else 'yes'}")
    print(f"saved output: {output_path}")


if __name__ == "__main__":
    main()
