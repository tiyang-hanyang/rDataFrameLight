import argparse
from pathlib import Path

import awkward as ak
import numpy as np
import uproot

from signal_class_config import DEFAULT_SCHEMA_NAME, get_signal_schema, schema_help_text


DEFAULT_INT = -1
DEFAULT_FLOAT = -999.0
DEFAULT_STRING = ""


def parse_args():
    parser = argparse.ArgumentParser(
        description="Append signal-classification prediction branches from a structured .npy file to a new ROOT file."
    )
    parser.add_argument("--input", required=True, help="Input ROOT file.")
    parser.add_argument("--prediction", required=True, help="Structured .npy produced by predictSignalRootToNpy.py.")
    parser.add_argument("--output", required=True, help="Output ROOT file.")
    parser.add_argument("--tree-name", default="Events", help="Input/output tree name.")
    parser.add_argument("--prefix", default="Signal", help="Prefix for appended branches.")
    parser.add_argument(
        "--class-schema",
        default=None,
        help="Optional legacy class schema. If omitted, score fields are inferred from the prediction file.",
    )
    parser.add_argument("--step-size", default="100 MB", help="Chunk size passed to uproot.iterate.")
    parser.add_argument(
        "--skip-copy-objects",
        action="store_true",
        help="Do not copy non-tree top-level ROOT objects such as genWeightSum.",
    )
    return parser.parse_args()


def infer_score_fields(prediction):
    known_fields = {"run", "luminosityBlock", "event", "prediction", "best_score", "predicted_class_name"}
    return [field for field in (prediction.dtype.names or ()) if field not in known_fields]


def resolve_tree_name(input_path, requested_tree_name):
    with uproot.open(input_path) as root_file:
        if requested_tree_name in root_file:
            return requested_tree_name

        tree_candidates = []
        for key in root_file.keys(cycle=False):
            try:
                obj = root_file[key]
            except Exception:
                continue
            if hasattr(obj, "num_entries") and hasattr(obj, "arrays"):
                tree_candidates.append(key)

    if len(tree_candidates) == 1:
        return tree_candidates[0]
    if not tree_candidates:
        raise RuntimeError(
            f"No TTree-like object found in {input_path}. Available keys do not include '{requested_tree_name}'."
        )
    raise RuntimeError(
        f"Requested tree '{requested_tree_name}' not found in {input_path}. "
        f"Available TTrees: {', '.join(tree_candidates)}"
    )


def load_prediction_table(prediction_path, score_fields):
    prediction = np.load(prediction_path, allow_pickle=False)
    required_fields = {"run", "luminosityBlock", "event", "prediction", *score_fields}
    missing = required_fields.difference(prediction.dtype.names or ())
    if missing:
        raise RuntimeError("Prediction file is missing required fields: " + ", ".join(sorted(missing)))

    table = {}
    for row in prediction:
        key = (int(row["run"]), int(row["luminosityBlock"]), int(row["event"]))
        table[key] = row
    return table, prediction


def build_new_columns(arrays, table, prefix, score_fields):
    runs = ak.to_numpy(arrays["run"])
    lumis = ak.to_numpy(arrays["luminosityBlock"])
    events = ak.to_numpy(arrays["event"])
    n_events = len(runs)

    prediction = np.full(n_events, DEFAULT_INT, dtype=np.int32)
    best_score = np.full(n_events, DEFAULT_FLOAT, dtype=np.float32)
    max_name_length = max([len(DEFAULT_STRING)] + [len(branch) for branch in score_fields])
    predicted_class_name = np.full(n_events, DEFAULT_STRING, dtype=f"<U{max_name_length}")
    scores = {branch: np.full(n_events, DEFAULT_FLOAT, dtype=np.float32) for branch in score_fields}

    matched = 0
    for index, key in enumerate(zip(runs, lumis, events)):
        row = table.get((int(key[0]), int(key[1]), int(key[2])))
        if row is None:
            continue
        matched += 1
        prediction[index] = int(row["prediction"])
        if "best_score" in row.dtype.names:
            best_score[index] = np.float32(row["best_score"])
        elif score_fields:
            best_score[index] = np.float32(max(float(row[branch]) for branch in score_fields))
        if "predicted_class_name" in row.dtype.names:
            predicted_class_name[index] = str(row["predicted_class_name"])
        for branch in score_fields:
            scores[branch][index] = np.float32(row[branch])

    new_columns = {
        f"{prefix}_prediction": prediction,
        f"{prefix}_bestScore": best_score,
        f"{prefix}_predictedClassName": ak.Array(predicted_class_name.tolist()),
        f"{prefix}_isMatched": prediction >= 0,
    }
    for branch in score_fields:
        new_columns[f"{prefix}_score_{branch}"] = scores[branch]
    return new_columns, matched


def main():
    args = parse_args()

    input_path = Path(args.input)
    prediction_path = Path(args.prediction)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    resolved_tree_name = resolve_tree_name(input_path, args.tree_name)
    score_fields = None

    table, prediction = load_prediction_table(prediction_path, [])
    score_fields = infer_score_fields(prediction)
    if args.class_schema is not None:
        expected_score_fields = list(get_signal_schema(args.class_schema).class_names)
        if expected_score_fields != score_fields:
            raise RuntimeError(
                f"Prediction file score fields are {score_fields}, but schema '{args.class_schema}' expects "
                f"{expected_score_fields}."
            )
    total_input_events = 0
    total_matched_events = 0
    has_written = False

    with uproot.open(input_path) as input_file, uproot.recreate(output_path) as output_file:
        if not args.skip_copy_objects:
            for key in input_file.keys(cycle=False):
                if key == resolved_tree_name:
                    continue
                output_file[key] = input_file[key]

        for arrays in uproot.iterate(f"{input_path}:{resolved_tree_name}", library="ak", step_size=args.step_size):
            total_input_events += len(arrays["run"])
            new_columns, matched = build_new_columns(arrays, table, args.prefix, score_fields)
            total_matched_events += matched

            merged = {field: arrays[field] for field in arrays.fields}
            merged.update(new_columns)

            if not has_written:
                output_file[resolved_tree_name] = merged
                has_written = True
            else:
                output_file[resolved_tree_name].extend(merged)

    print(f"tree name: {resolved_tree_name}")
    print(f"prediction rows loaded: {len(table)}")
    print(f"input events written: {total_input_events}")
    print(f"matched events filled: {total_matched_events}")
    print(f"copied extra ROOT objects: {'no' if args.skip_copy_objects else 'yes'}")
    print(f"saved output: {output_path}")


if __name__ == "__main__":
    main()
