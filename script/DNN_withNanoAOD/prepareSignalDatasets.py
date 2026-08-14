import argparse
import json
from pathlib import Path
import warnings

import numpy as np

from signalExtractionInputFromNanoAOD import ALL_FIELD_NAMES as EXTRACTED_FIELD_NAMES, FEATURE_FIELD_NAMES
from signal_class_config import DEFAULT_SCHEMA_NAME, get_signal_schema, schema_help_text

FEATURE_DIM = len(FEATURE_FIELD_NAMES)
EXTRACTED_WEIGHT_INDEX = FEATURE_DIM
EXTRACTED_SIGN_INDEX = FEATURE_DIM + 1
WEIGHT_INDEX = FEATURE_DIM
LABEL_INDEX = FEATURE_DIM + 1
SIGN_INDEX = FEATURE_DIM + 2
CHANNEL_INDEX = FEATURE_DIM + 3
FINAL_FIELD_NAMES = FEATURE_FIELD_NAMES + ["weight", "label", "isSameSign"]
MIXED_FIELD_NAMES = FINAL_FIELD_NAMES + ["channelId"]
DEFAULT_MIX_REPORT_NAME = "mix_report.json"


def structured_to_matrix(array):
    if array.dtype.names is None:
        return array.astype(np.float32)
    dtype_names = array.dtype.names or ()
    preferred_field_names = MIXED_FIELD_NAMES if "label" in dtype_names else EXTRACTED_FIELD_NAMES
    field_names = [field_name for field_name in preferred_field_names if field_name in array.dtype.names]
    missing = [field_name for field_name in preferred_field_names if field_name not in dtype_names]
    if missing:
        raise RuntimeError(
            "Dataset is missing required fields: " + ", ".join(missing)
        )
    columns = [array[field_name].astype(np.float32) for field_name in field_names]
    return np.column_stack(columns).astype(np.float32)


def matrix_to_structured(array):
    dtype = [(field_name, np.float32) for field_name in MIXED_FIELD_NAMES]
    structured = np.empty(array.shape[0], dtype=dtype)
    if array.size == 0:
        return structured
    for field_index, field_name in enumerate(MIXED_FIELD_NAMES):
        structured[field_name] = array[:, field_index]
    return structured


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Merge per-sample signal datasets and split OS into train/val/test while keeping SS as SRtest. "
            "This is the main dataset-preparation entry point."
        )
    )
    parser.add_argument(
        "--inputs",
        nargs="+",
        required=True,
        help="Input .npy datasets produced by signalExtractionInputFromNanoAOD.py",
    )
    parser.add_argument("--output-dir", required=True, help="Output directory.")
    parser.add_argument("--val-frac", type=float, default=0.15, help="Validation fraction within OS.")
    parser.add_argument("--test-frac", type=float, default=0.15, help="Test fraction within OS.")
    parser.add_argument("--seed", type=int, default=12345, help="Random seed.")
    parser.add_argument(
        "--weighting-mode",
        choices=["balanced", "raw"],
        default="balanced",
        help=(
            "'balanced' is the main workflow: rebalance OS_train/OS_val weights by class and subgroup. "
            "'raw' only splits and keeps original weights."
        ),
    )
    parser.add_argument(
        "--class-schema",
        default=DEFAULT_SCHEMA_NAME,
        help=schema_help_text(),
    )
    parser.add_argument(
        "--save-debug-files",
        action="store_true",
        help="Also save optional debug files such as raw OS_train/OS_val copies and balance reports.",
    )
    return parser.parse_args()


def load_and_concat(paths):
    channel_map = {Path(path).stem: index for index, path in enumerate(paths)}
    arrays = []
    input_records = []
    for path in paths:
        raw_array = np.load(path, allow_pickle=False)
        matrix = structured_to_matrix(raw_array)
        channel_name = Path(path).stem
        channel_id = channel_map[channel_name]
        input_records.append(
            {
                "path": str(Path(path).resolve()),
                "channel_name": channel_name,
                "channel_id": int(channel_id),
                "input_fields": list(raw_array.dtype.names or EXTRACTED_FIELD_NAMES),
                "n_events": int(matrix.shape[0]),
                "weight_sum": float(np.sum(matrix[:, EXTRACTED_WEIGHT_INDEX], dtype=np.float64)) if matrix.size > 0 else 0.0,
                "abs_weight_sum": float(np.sum(np.abs(matrix[:, EXTRACTED_WEIGHT_INDEX]), dtype=np.float64)) if matrix.size > 0 else 0.0,
                "same_sign_events": int(np.sum(matrix[:, EXTRACTED_SIGN_INDEX] > 0.5)) if matrix.size > 0 else 0,
                "opposite_sign_events": int(np.sum(matrix[:, EXTRACTED_SIGN_INDEX] <= 0.5)) if matrix.size > 0 else 0,
            }
        )
        if matrix.size == 0:
            continue
        weights = matrix[:, EXTRACTED_WEIGHT_INDEX : EXTRACTED_WEIGHT_INDEX + 1]
        signs = matrix[:, EXTRACTED_SIGN_INDEX : EXTRACTED_SIGN_INDEX + 1]
        labels = np.full((matrix.shape[0], 1), -1.0, dtype=np.float32)
        channel_column = np.full((matrix.shape[0], 1), channel_id, dtype=np.float32)
        features = matrix[:, :FEATURE_DIM]
        arrays.append(np.concatenate([features, weights, labels, signs, channel_column], axis=1).astype(np.float32))
    non_empty = [array for array in arrays if array.size > 0]
    if not non_empty:
        return np.empty((0, CHANNEL_INDEX + 1), dtype=np.float32), channel_map, input_records
    return np.concatenate(non_empty, axis=0).astype(np.float32), channel_map, input_records


def resolve_input_paths(inputs):
    resolved = []
    seen = set()
    for raw_input in inputs:
        path = Path(raw_input)
        if path.is_dir():
            candidates = sorted(candidate for candidate in path.rglob("*.npy") if candidate.is_file())
        elif path.is_file():
            candidates = [path]
        else:
            raise RuntimeError(f"Input path does not exist: {raw_input}")

        for candidate in candidates:
            candidate_resolved = candidate.resolve()
            candidate_key = str(candidate_resolved)
            if candidate_key in seen:
                continue
            seen.add(candidate_key)
            resolved.append(candidate_resolved)

    if not resolved:
        raise RuntimeError("No input .npy files were found from --inputs.")
    return resolved


def stratified_split(dataset, val_frac, test_frac, seed):
    rng = np.random.default_rng(seed)
    labels = dataset[:, LABEL_INDEX].astype(int)
    channels = dataset[:, CHANNEL_INDEX].astype(int)
    group_keys = np.unique(np.column_stack([labels, channels]), axis=0)

    train_parts = []
    val_parts = []
    test_parts = []

    for label, channel in group_keys:
        group_rows = dataset[(labels == label) & (channels == channel)]
        if group_rows.size == 0:
            continue
        permuted = group_rows[rng.permutation(group_rows.shape[0])]
        n_total = permuted.shape[0]
        n_val = int(round(n_total * val_frac))
        n_test = int(round(n_total * test_frac))
        if n_val + n_test >= n_total and n_total > 1:
            overflow = n_val + n_test - (n_total - 1)
            n_test = max(0, n_test - overflow)
        n_train = max(0, n_total - n_val - n_test)

        train_parts.append(permuted[:n_train])
        val_parts.append(permuted[n_train : n_train + n_val])
        test_parts.append(permuted[n_train + n_val : n_train + n_val + n_test])

    def merge(parts):
        valid_parts = [part for part in parts if part.size > 0]
        if not valid_parts:
            return np.empty((0, dataset.shape[1]), dtype=np.float32)
        merged = np.concatenate(valid_parts, axis=0).astype(np.float32)
        return merged[rng.permutation(merged.shape[0])]

    return merge(train_parts), merge(val_parts), merge(test_parts)


def save_dataset(path, array):
    path.parent.mkdir(parents=True, exist_ok=True)
    np.save(path, matrix_to_structured(array.astype(np.float32)))


def save_channel_map(path, channel_map):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(channel_map, handle, indent=2, sort_keys=True)


def save_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)


def balance_group_name(channel_name, class_name):
    if channel_name in {"ttbarDL", "ttbarSL", "ttbar"}:
        return "ttbar"
    if channel_name in {"TTBB_DL", "TTBB_SL", "TTBB"}:
        return "TTBB"
    if channel_name in {"TTZ_low", "TTZ_high", "TTZ"}:
        return "TTZ"
    return channel_name if class_name not in {"ttbar_ttbb", "tt_b"} else channel_name


def canonical_channel_map(channel_map, schema):
    canonical = {}
    skipped = {}
    for original_name, channel_id in channel_map.items():
        try:
            channel_name, class_name = schema.identify_channel(original_name)
        except RuntimeError as exc:
            skipped[int(channel_id)] = {
                "original_name": original_name,
                "reason": str(exc),
            }
            continue
        canonical[int(channel_id)] = {
            "original_name": original_name,
            "channel_name": channel_name,
            "class_name": class_name,
            "label_index": int(schema.class_to_index[class_name]),
            "balance_group": (
                next(
                    (
                        channel.balance_group
                        for channel in schema.channels
                        if channel.channel_name == channel_name and channel.class_name == class_name and channel.balance_group
                    ),
                    None,
                )
                or balance_group_name(channel_name, class_name)
            ),
        }
    return canonical, skipped


def relabel_dataset_from_channels(dataset, channel_info):
    if dataset.size == 0:
        return dataset.astype(np.float32)

    relabeled = dataset.astype(np.float32).copy()
    channel_ids = relabeled[:, CHANNEL_INDEX].astype(np.int64)
    relabeled_labels = np.empty(relabeled.shape[0], dtype=np.float32)

    for row_index, channel_id in enumerate(channel_ids):
        info = channel_info.get(int(channel_id))
        if info is None:
            raise RuntimeError(f"Missing channel metadata for channelId={channel_id}.")
        relabeled_labels[row_index] = np.float32(info["label_index"])

    relabeled[:, LABEL_INDEX] = relabeled_labels
    return relabeled


def filter_dataset_to_known_channels(dataset, channel_info):
    if dataset.size == 0:
        return dataset.astype(np.float32)
    known_channel_ids = set(int(channel_id) for channel_id in channel_info)
    channel_ids = dataset[:, CHANNEL_INDEX].astype(np.int64)
    keep_mask = np.array([int(channel_id) in known_channel_ids for channel_id in channel_ids], dtype=bool)
    return dataset[keep_mask].astype(np.float32)


def filter_dataset_by_schema_split(dataset, channel_info, schema, split_name):
    if dataset.size == 0:
        return dataset.astype(np.float32)

    channel_ids = dataset[:, CHANNEL_INDEX].astype(np.int64)
    keep_mask = np.zeros(dataset.shape[0], dtype=bool)
    for row_index, channel_id in enumerate(channel_ids):
        info = channel_info.get(int(channel_id))
        if info is None:
            continue
        channel_name = info["channel_name"]
        if split_name == "os":
            keep_mask[row_index] = schema.is_allowed_in_os(channel_name)
        elif split_name == "sr":
            keep_mask[row_index] = schema.is_allowed_in_sr(channel_name)
        else:
            raise RuntimeError(f"Unknown split_name={split_name}")
    return dataset[keep_mask].astype(np.float32)


def compute_balanced_training_weights(dataset, channel_info):
    if dataset.size == 0:
        return np.empty((0,), dtype=np.float32), {}

    labels = dataset[:, LABEL_INDEX].astype(np.int64)
    channel_ids = dataset[:, CHANNEL_INDEX].astype(np.int64)
    raw_weights = np.abs(dataset[:, WEIGHT_INDEX].astype(np.float64))

    class_groups = {}
    for row_index, (label, channel_id) in enumerate(zip(labels, channel_ids)):
        info = channel_info.get(int(channel_id))
        if info is None:
            continue
        group_key = (int(label), info["balance_group"])
        class_groups.setdefault(group_key, []).append(row_index)

    group_sums = {
        group_key: float(np.sum(raw_weights[row_indices], dtype=np.float64))
        for group_key, row_indices in class_groups.items()
    }

    class_targets = {}
    for label in np.unique(labels):
        class_group_keys = [key for key in group_sums if key[0] == int(label) and group_sums[key] > 0.0]
        if not class_group_keys:
            continue
        class_targets[int(label)] = 1.0 / float(len(class_group_keys))

    balanced = raw_weights.copy()
    scale_report = {}
    for group_key, row_indices in class_groups.items():
        group_sum = group_sums[group_key]
        if group_sum <= 0.0 or group_key[0] not in class_targets:
            scale = 0.0
        else:
            scale = class_targets[group_key[0]] / group_sum
        balanced[row_indices] *= scale
        scale_report[f"{group_key[0]}::{group_key[1]}"] = {
            "raw_sum": group_sum,
            "scale": scale,
            "balanced_sum": float(np.sum(balanced[row_indices], dtype=np.float64)),
            "n_events": int(len(row_indices)),
        }

    return balanced.astype(np.float32), scale_report


def replace_weight_column(dataset, new_weights):
    updated = dataset.astype(np.float32).copy()
    updated[:, WEIGHT_INDEX] = new_weights.astype(np.float32)
    return updated


def normalize_channel_info_for_json(channel_info):
    return {str(channel_id): info for channel_id, info in sorted(channel_info.items())}


def count_values(values):
    if values.size == 0:
        return {}
    unique, counts = np.unique(values, return_counts=True)
    return {int(value): int(count) for value, count in zip(unique, counts)}


def sum_weights_by_key(keys, weights):
    summary = {}
    for key in np.unique(keys):
        mask = keys == key
        summary[int(key)] = float(np.sum(weights[mask], dtype=np.float64))
    return summary


def summarize_dataset_for_report(dataset, channel_info, label_names):
    matrix = structured_to_matrix(dataset)
    summary = {
        "n_events": int(matrix.shape[0]),
        "feature_dim": FEATURE_DIM,
        "feature_names": list(FEATURE_FIELD_NAMES),
        "mixed_field_names": list(MIXED_FIELD_NAMES),
        "label_names": list(label_names),
        "weight_sum": 0.0,
        "abs_weight_sum": 0.0,
        "same_sign_events": 0,
        "opposite_sign_events": 0,
        "label_event_counts": {},
        "label_weight_sums": {},
        "channel_event_counts": {},
        "channel_weight_sums": {},
    }
    if matrix.size == 0:
        return summary

    weights = matrix[:, WEIGHT_INDEX].astype(np.float64)
    labels = matrix[:, LABEL_INDEX].astype(np.int64)
    signs = matrix[:, SIGN_INDEX].astype(np.int64)
    channels = matrix[:, CHANNEL_INDEX].astype(np.int64)
    summary["weight_sum"] = float(np.sum(weights, dtype=np.float64))
    summary["abs_weight_sum"] = float(np.sum(np.abs(weights), dtype=np.float64))
    summary["same_sign_events"] = int(np.sum(signs > 0))
    summary["opposite_sign_events"] = int(np.sum(signs <= 0))

    label_counts = count_values(labels)
    label_weight_sums = sum_weights_by_key(labels, weights)
    summary["label_event_counts"] = {
        label_names[label_index]: label_counts.get(label_index, 0)
        for label_index in range(len(label_names))
    }
    summary["label_weight_sums"] = {
        label_names[label_index]: label_weight_sums.get(label_index, 0.0)
        for label_index in range(len(label_names))
    }

    channel_counts = count_values(channels)
    channel_weight_sums = sum_weights_by_key(channels, weights)
    for channel_id, info in sorted(channel_info.items()):
        channel_name = info["original_name"]
        summary["channel_event_counts"][channel_name] = channel_counts.get(channel_id, 0)
        summary["channel_weight_sums"][channel_name] = channel_weight_sums.get(channel_id, 0.0)
    return summary


def enrich_input_records(input_records, channel_info, label_names):
    enriched = []
    for record in input_records:
        info = channel_info.get(int(record["channel_id"]))
        if info is None:
            enriched.append(
                {
                    **record,
                    "included_in_schema": False,
                }
            )
            continue
        enriched.append(
            {
                **record,
                "included_in_schema": True,
                "canonical_channel_name": info["channel_name"],
                "class_name": info["class_name"],
                "label_index": info["label_index"],
                "label_name": label_names[info["label_index"]],
                "balance_group": info["balance_group"],
            }
        )
    return enriched


def extract_common_input_fields(input_records):
    if not input_records:
        return list(EXTRACTED_FIELD_NAMES), []

    common_fields = input_records[0].get("input_fields", [])
    compact_records = []
    for record in input_records:
        compact_record = dict(record)
        record_fields = compact_record.get("input_fields", [])
        if record_fields == common_fields:
            compact_record.pop("input_fields", None)
        compact_records.append(compact_record)
    return common_fields, compact_records


def build_weighting_report(before_dataset, after_dataset, channel_info, scale_report, label_names):
    return {
        "before": summarize_dataset_for_report(before_dataset, channel_info, label_names),
        "after": summarize_dataset_for_report(after_dataset, channel_info, label_names),
        "group_scales": scale_report,
    }


def build_mix_report(
    output_dir,
    schema_name,
    label_names,
    weighting_mode,
    schema,
    channel_map,
    channel_info,
    skipped_channels,
    input_records,
    train_set,
    val_set,
    test_set,
    same_sign,
    train_set_raw=None,
    val_set_raw=None,
    balance_report=None,
):
    common_input_fields, compact_input_records = extract_common_input_fields(input_records)
    report = {
        "output_dir": str(output_dir.resolve()),
        "label_names": list(label_names),
        "feature_names": list(FEATURE_FIELD_NAMES),
        "extracted_field_names": list(EXTRACTED_FIELD_NAMES),
        "mixed_field_names": list(MIXED_FIELD_NAMES),
        "common_input_fields": common_input_fields,
        "labeling": {
            "policy": (
                "Extraction remains classification-agnostic. "
                "Labels are assigned only in prepareSignalDatasets.py from channelId during mix/split."
            ),
            "label_config_name": schema_name,
            "allowed_os_channels": None if schema.allowed_os_channels is None else list(schema.allowed_os_channels),
            "allowed_sr_channels": None if schema.allowed_sr_channels is None else list(schema.allowed_sr_channels),
            "channel_map": channel_map,
            "channel_info": normalize_channel_info_for_json(channel_info),
            "skipped_channels": {str(channel_id): info for channel_id, info in sorted(skipped_channels.items())},
        },
        "inputs": compact_input_records,
        "splits": {
            "OS_train": summarize_dataset_for_report(train_set, channel_info, label_names),
            "OS_val": summarize_dataset_for_report(val_set, channel_info, label_names),
            "OS_test": summarize_dataset_for_report(test_set, channel_info, label_names),
            "SRtest": summarize_dataset_for_report(same_sign, channel_info, label_names),
        },
        "reweighting": {
            "mode": weighting_mode,
            "train": None,
            "val": None,
        },
    }
    if weighting_mode == "balanced" and train_set_raw is not None and val_set_raw is not None and balance_report is not None:
        report["reweighting"]["train"] = build_weighting_report(
            train_set_raw,
            train_set,
            channel_info,
            balance_report["train"],
            label_names,
        )
        report["reweighting"]["val"] = build_weighting_report(
            val_set_raw,
            val_set,
            channel_info,
            balance_report["val"],
            label_names,
        )
    return report


def load_json(path):
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def resolve_dataset_report_path(dataset_path, explicit_report_path=None):
    if explicit_report_path is not None:
        return Path(explicit_report_path)
    candidate = Path(dataset_path).resolve().parent / DEFAULT_MIX_REPORT_NAME
    if candidate.exists():
        return candidate
    raise RuntimeError(
        f"Unable to find {DEFAULT_MIX_REPORT_NAME} next to dataset {dataset_path}. "
        "Pass --dataset-report explicitly."
    )


def validate_dataset_report(report, context):
    if report.get("feature_names") != list(FEATURE_FIELD_NAMES):
        raise RuntimeError(f"{context}: feature_names do not match the code-side FEATURE_FIELD_NAMES contract.")
    if report.get("mixed_field_names") != list(MIXED_FIELD_NAMES):
        raise RuntimeError(f"{context}: mixed_field_names do not match the expected mixed dataset fields.")
    label_names = report.get("label_names")
    if not isinstance(label_names, list) or not label_names:
        raise RuntimeError(f"{context}: report is missing a non-empty label_names list.")
    if "labeling" not in report or "channel_info" not in report["labeling"]:
        raise RuntimeError(f"{context}: report is missing labeling.channel_info.")


def validate_mixed_dataset_file(dataset_path, report, expected_sign=None):
    dataset_path = Path(dataset_path).resolve()
    validate_dataset_report(report, str(dataset_path))
    raw_dataset = np.load(dataset_path, allow_pickle=False)
    actual_fields = list(raw_dataset.dtype.names or [])
    if actual_fields != list(MIXED_FIELD_NAMES):
        raise RuntimeError(
            f"{dataset_path}: field mismatch. Expected {list(MIXED_FIELD_NAMES)}, got {actual_fields}."
        )

    matrix = structured_to_matrix(raw_dataset)
    if matrix.shape[1] != len(MIXED_FIELD_NAMES):
        raise RuntimeError(
            f"{dataset_path}: expected {len(MIXED_FIELD_NAMES)} columns, got {matrix.shape[1]}."
        )

    if matrix.size == 0:
        return {"path": str(dataset_path), "n_events": 0, "status": "ok"}

    labels = matrix[:, LABEL_INDEX].astype(np.int64)
    channel_ids = matrix[:, CHANNEL_INDEX].astype(np.int64)
    signs = matrix[:, SIGN_INDEX].astype(np.int64)
    label_names = report["label_names"]
    if np.any(labels < 0) or np.any(labels >= len(label_names)):
        raise RuntimeError(f"{dataset_path}: found labels outside [0, {len(label_names) - 1}].")

    channel_info = {int(key): value for key, value in report["labeling"]["channel_info"].items()}
    missing_channels = sorted(set(channel_ids.tolist()) - set(channel_info))
    if missing_channels:
        raise RuntimeError(f"{dataset_path}: missing channel metadata for channelIds {missing_channels}.")

    expected_labels = np.array([channel_info[int(channel_id)]["label_index"] for channel_id in channel_ids], dtype=np.int64)
    mismatch = np.flatnonzero(labels != expected_labels)
    if mismatch.size > 0:
        first = int(mismatch[0])
        raise RuntimeError(
            f"{dataset_path}: label/channel mismatch at row {first}. "
            f"label={int(labels[first])}, expected={int(expected_labels[first])}, channelId={int(channel_ids[first])}."
        )

    if expected_sign == "os" and np.any(signs > 0):
        raise RuntimeError(f"{dataset_path}: expected opposite-sign dataset but found same-sign events.")
    if expected_sign == "ss" and np.any(signs <= 0):
        raise RuntimeError(f"{dataset_path}: expected same-sign dataset but found opposite-sign events.")

    return {
        "path": str(dataset_path),
        "n_events": int(matrix.shape[0]),
        "label_counts": count_values(labels),
        "channel_counts": count_values(channel_ids),
        "status": "ok",
    }


def summarize(name, dataset):
    dataset = structured_to_matrix(dataset)
    labels = dataset[:, LABEL_INDEX].astype(int) if dataset.size > 0 else np.array([], dtype=int)
    unique, counts = np.unique(labels, return_counts=True)
    count_map = {int(label): int(count) for label, count in zip(unique, counts)}
    channels = dataset[:, CHANNEL_INDEX].astype(int) if dataset.size > 0 else np.array([], dtype=int)
    unique_channels, channel_counts = np.unique(channels, return_counts=True)
    channel_map = {int(channel): int(count) for channel, count in zip(unique_channels, channel_counts)}
    print(f"{name}: {dataset.shape[0]} events, label counts = {count_map}, channel counts = {channel_map}")


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    input_paths = resolve_input_paths(args.inputs)
    schema = get_signal_schema(args.class_schema)
    label_names = list(schema.class_names)

    dataset, channel_map, input_records = load_and_concat(input_paths)
    channel_info, skipped_channels = canonical_channel_map(channel_map, schema)
    if skipped_channels:
        skipped_names = ", ".join(sorted(info["original_name"] for info in skipped_channels.values()))
        warnings.warn(
            f"Skipping inputs not used by schema '{args.class_schema}': {skipped_names}",
            stacklevel=2,
        )
    enriched_inputs = enrich_input_records(input_records, channel_info, label_names)
    dataset = filter_dataset_to_known_channels(dataset, channel_info)
    dataset = relabel_dataset_from_channels(dataset, channel_info)
    same_sign = filter_dataset_by_schema_split(dataset[dataset[:, SIGN_INDEX] > 0.5], channel_info, schema, "sr")
    opposite_sign = filter_dataset_by_schema_split(dataset[dataset[:, SIGN_INDEX] <= 0.5], channel_info, schema, "os")

    train_set, val_set, test_set = stratified_split(
        opposite_sign,
        val_frac=args.val_frac,
        test_frac=args.test_frac,
        seed=args.seed,
    )
    train_set_raw = train_set.copy()
    val_set_raw = val_set.copy()

    balance_report = None
    if args.weighting_mode == "balanced":
        train_weights, train_report = compute_balanced_training_weights(train_set_raw, channel_info)
        val_weights, val_report = compute_balanced_training_weights(val_set_raw, channel_info)
        train_set = replace_weight_column(train_set, train_weights)
        val_set = replace_weight_column(val_set, val_weights)
        balance_report = {
            "train": train_report,
            "val": val_report,
        }

    save_dataset(output_dir / "OS_train.npy", train_set)
    save_dataset(output_dir / "OS_val.npy", val_set)
    save_dataset(output_dir / "OS_test.npy", test_set)
    save_dataset(output_dir / "SRtest.npy", same_sign)
    save_channel_map(output_dir / "channel_map.json", channel_map)
    save_json(
        output_dir / DEFAULT_MIX_REPORT_NAME,
        build_mix_report(
            output_dir=output_dir,
            schema_name=args.class_schema,
            label_names=label_names,
            weighting_mode=args.weighting_mode,
            schema=schema,
            channel_map=channel_map,
            channel_info=channel_info,
            skipped_channels=skipped_channels,
            input_records=enriched_inputs,
            train_set=train_set,
            val_set=val_set,
            test_set=test_set,
            same_sign=same_sign,
            train_set_raw=train_set_raw,
            val_set_raw=val_set_raw,
            balance_report=balance_report,
        ),
    )

    if args.save_debug_files and args.weighting_mode == "balanced":
        save_dataset(output_dir / "OS_train_raw.npy", train_set_raw)
        save_dataset(output_dir / "OS_val_raw.npy", val_set_raw)
        save_json(output_dir / "channel_balance_groups.json", channel_info)
        save_json(output_dir / "balance_weight_report.json", balance_report)

    summarize(f"OS_train ({args.weighting_mode})", matrix_to_structured(train_set))
    summarize(f"OS_val ({args.weighting_mode})", matrix_to_structured(val_set))
    summarize("OS_test", matrix_to_structured(test_set))
    summarize("SRtest", matrix_to_structured(same_sign))
    print(f"resolved input files: {len(input_paths)}")
    print(f"saved mix report: {output_dir / DEFAULT_MIX_REPORT_NAME}")
    print(f"saved output dir: {output_dir}")


if __name__ == "__main__":
    main()
