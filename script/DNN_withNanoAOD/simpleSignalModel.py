import argparse
from pathlib import Path
import warnings

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset

from matplotlib import pyplot as plt

from feature_grouping import (
    apply_feature_mask,
    get_active_feature_names,
    get_group_names,
    resolve_feature_groups,
    select_active_features,
)
from prepareSignalDatasets import (
    CHANNEL_INDEX,
    DEFAULT_MIX_REPORT_NAME,
    FEATURE_DIM,
    LABEL_INDEX,
    MIXED_FIELD_NAMES,
    SIGN_INDEX,
    WEIGHT_INDEX,
    load_json,
    resolve_dataset_report_path,
    structured_to_matrix,
    validate_dataset_report,
    validate_mixed_dataset_file,
)
from signalExtractionInputFromNanoAOD import FEATURE_FIELD_NAMES
from signal_class_config import DEFAULT_SCHEMA_NAME, get_signal_schema
from training_utils import plot_confusion_matrix_fractional, plot_loss_curves


CLASS_NAMES = list(get_signal_schema().class_names)


class SimpleMLP(nn.Module):
    def __init__(self, input_dim=FEATURE_DIM, hidden_dims=None, dropout=0.2, n_classes=None):
        super().__init__()
        if hidden_dims is None:
            hidden_dims = [256, 128, 64]
        if n_classes is None:
            n_classes = len(CLASS_NAMES)

        layers = []
        in_dim = input_dim
        for hidden_dim in hidden_dims:
            layers.append(nn.Linear(in_dim, hidden_dim))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(dropout))
            in_dim = hidden_dim
        layers.append(nn.Linear(in_dim, n_classes))
        self.network = nn.Sequential(*layers)

    def forward(self, x):
        return self.network(x)


def extract_model_dropout_rate(model):
    dropout_layers = [layer.p for layer in model.modules() if isinstance(layer, nn.Dropout)]
    if not dropout_layers:
        return 0.0
    first = float(dropout_layers[0])
    if any(abs(float(rate) - first) > 1e-12 for rate in dropout_layers[1:]):
        return dropout_layers
    return first


def print_training_start_summary(args, model, class_names, report_path):
    print("Training configuration:")
    print(f"  dataset report: {report_path}")
    print(f"  train dataset: {args.train}")
    print(f"  val dataset: {args.val}")
    print(f"  output dir: {args.output_dir}")
    print(f"  epochs: {args.epochs}")
    print(f"  batch size: {args.batch_size}")
    print(f"  learning rate: {args.lr}")
    print(f"  training weight mode: {args.training_weight_mode}")
    print(f"  weight targets: {len(getattr(args, 'weight_target', []) or [])}")
    print(f"  n_classes: {len(class_names)}")
    print(f"  class names: {', '.join(class_names)}")
    if getattr(args, "label_subset", None) is not None:
        print(f"  label subset: {', '.join(args.label_subset)}")
    if getattr(args, "feature_groups", None) is not None:
        print(f"  feature groups: {', '.join(resolve_feature_groups(args.feature_groups))}")
    if getattr(args, "feature_layout", None) is not None:
        print(f"  feature layout: {args.feature_layout}")
    print(f"  dropout rate: {extract_model_dropout_rate(model)}")
    print("Model architecture:")
    print(model)


def loadDS(file_path):
    dataset = structured_to_matrix(np.load(file_path))
    x = dataset[:, :FEATURE_DIM].astype(np.float32)
    weights = dataset[:, WEIGHT_INDEX].astype(np.float32)
    labels = dataset[:, LABEL_INDEX].astype(np.int64)
    is_same_sign = dataset[:, SIGN_INDEX].astype(np.int64)
    channel_ids = (
        dataset[:, CHANNEL_INDEX].astype(np.int64)
        if dataset.shape[1] > CHANNEL_INDEX
        else np.full(dataset.shape[0], -1, dtype=np.int64)
    )
    return x, weights, labels, is_same_sign, channel_ids


def build_loader(x, y, weights, batch_size, shuffle):
    if weights is None:
        dataset = TensorDataset(
            torch.from_numpy(x).float(),
            torch.from_numpy(y).long(),
        )
    else:
        dataset = TensorDataset(
            torch.from_numpy(x).float(),
            torch.from_numpy(y).long(),
            torch.from_numpy(weights).float(),
        )
    return DataLoader(dataset, batch_size=batch_size, shuffle=shuffle, drop_last=False)


def normalize_event_weights_by_class(weights, labels):
    normalized = weights.astype(np.float32).copy()
    unique_labels = np.unique(labels)
    for label in unique_labels:
        class_mask = labels == label
        class_sum = np.sum(normalized[class_mask], dtype=np.float64)
        if class_sum > 0:
            normalized[class_mask] /= np.float32(class_sum)
    return normalized


def make_non_negative_training_weights(weights):
    return np.abs(weights.astype(np.float32))


def build_balanced_training_weights(weights, labels):
    non_negative = make_non_negative_training_weights(weights)
    return normalize_event_weights_by_class(non_negative, labels)


def resolve_training_weights(weights, labels, mode):
    if mode == "balanced":
        return build_balanced_training_weights(weights, labels)
    if mode == "input":
        return weights.astype(np.float32)
    if mode == "abs":
        return make_non_negative_training_weights(weights)
    raise RuntimeError(f"Unknown training weight mode: {mode}")


def load_channel_info_from_report(dataset_report):
    raw_channel_info = dataset_report.get("labeling", {}).get("channel_info", {})
    return {int(key): value for key, value in raw_channel_info.items()}


def resolve_weight_targets(args):
    class_targets = {}
    channel_share_targets = {}
    for entry in getattr(args, "weight_target", []) or []:
        if "=" not in entry:
            raise RuntimeError(
                f"Invalid --weight-target entry '{entry}'. Use 'class=value' or 'class/channel=share'."
            )
        selector, value_text = entry.split("=", 1)
        try:
            value = float(value_text)
        except ValueError as exc:
            raise RuntimeError(f"Invalid numeric value '{value_text}' in --weight-target '{entry}'.") from exc
        if "/" in selector:
            class_name, channel_name = selector.split("/", 1)
            if not (0.0 <= value <= 1.0):
                raise RuntimeError(
                    f"Channel-share target '{entry}' must be in [0, 1] because it represents an in-class fraction."
                )
            channel_share_targets.setdefault(class_name, {})[channel_name] = value
        else:
            if value < 0.0:
                raise RuntimeError(f"Class target '{entry}' must be non-negative.")
            class_targets[selector] = value
    return class_targets, channel_share_targets


def apply_class_weight_targets(weights, labels, class_names, class_targets):
    if not class_targets:
        return weights.astype(np.float32).copy(), []

    effective = weights.astype(np.float32).copy()
    total_weight = float(np.sum(effective, dtype=np.float64))
    if total_weight <= 0.0:
        return effective, []

    unknown_classes = sorted(set(class_targets) - set(class_names))
    if unknown_classes:
        raise RuntimeError(f"Unknown class names in --weight-target: {', '.join(unknown_classes)}")

    class_scales = []
    present_labels = [int(label) for label in np.unique(labels)]
    target_sum_denominator = 0.0
    current_sums = {}
    for label in present_labels:
        class_mask = labels == label
        current_sum = float(np.sum(effective[class_mask], dtype=np.float64))
        current_sums[label] = current_sum
        if current_sum <= 0.0:
            continue
        target_sum_denominator += float(class_targets.get(class_names[label], 1.0))

    if target_sum_denominator <= 0.0:
        raise RuntimeError("All requested class targets are zero for the present classes.")

    for label in present_labels:
        class_name = class_names[label]
        current_sum = current_sums[label]
        if current_sum <= 0.0:
            continue
        target_relative = float(class_targets.get(class_name, 1.0))
        target_sum = total_weight * target_relative / target_sum_denominator
        scale = target_sum / current_sum
        effective[labels == label] *= np.float32(scale)
        class_scales.append(
            {
                "class_name": class_name,
                "before_sum": current_sum,
                "after_sum": float(np.sum(effective[labels == label], dtype=np.float64)),
                "target_relative": target_relative,
                "effective_scale": scale,
            }
        )
    return effective, class_scales


def apply_class_channel_share_targets(weights, labels, channel_ids, channel_info, class_names, channel_share_targets):
    if not channel_share_targets:
        return weights.astype(np.float32).copy(), []

    effective = weights.astype(np.float32).copy()
    summary = []
    channel_names = np.array(
        [str(channel_info.get(int(channel_id), {}).get("channel_name", f"channel_{int(channel_id)}")) for channel_id in channel_ids],
        dtype=object,
    )

    unknown_classes = sorted(set(channel_share_targets) - set(class_names))
    if unknown_classes:
        raise RuntimeError(f"Unknown class names in class/channel --weight-target: {', '.join(unknown_classes)}")

    for class_name, target_map in channel_share_targets.items():
        class_index = class_names.index(class_name)
        class_mask = labels == class_index
        if not np.any(class_mask):
            continue

        class_total = float(np.sum(effective[class_mask], dtype=np.float64))
        if class_total <= 0.0:
            continue

        class_channel_names = channel_names[class_mask]
        class_weights = effective[class_mask].astype(np.float64)
        unique_channels = sorted(set(class_channel_names.tolist()))
        current_channel_sums = {
            channel_name: float(np.sum(class_weights[class_channel_names == channel_name], dtype=np.float64))
            for channel_name in unique_channels
        }

        unknown_channels = sorted(set(target_map) - set(unique_channels))
        if unknown_channels:
            raise RuntimeError(
                f"Unknown channels for class '{class_name}' in --weight-target: {', '.join(unknown_channels)}"
            )

        fixed_share_sum = float(sum(target_map.values()))
        if fixed_share_sum > 1.0 + 1e-9:
            raise RuntimeError(
                f"Requested in-class shares for '{class_name}' sum to {fixed_share_sum:.6g}, larger than 1."
            )

        untargeted_channels = [channel_name for channel_name in unique_channels if channel_name not in target_map]
        remaining_share = max(0.0, 1.0 - fixed_share_sum)
        desired_channel_shares = dict(target_map)
        if untargeted_channels:
            untargeted_sum = float(sum(current_channel_sums[channel_name] for channel_name in untargeted_channels))
            if untargeted_sum <= 0.0 and remaining_share > 0.0:
                raise RuntimeError(
                    f"Cannot distribute remaining share inside '{class_name}' because untargeted channels have zero weight."
                )
            for channel_name in untargeted_channels:
                desired_channel_shares[channel_name] = (
                    remaining_share * current_channel_sums[channel_name] / untargeted_sum if untargeted_sum > 0.0 else 0.0
                )
        elif remaining_share > 1e-9:
            raise RuntimeError(
                f"Requested in-class shares for '{class_name}' sum to less than 1, but no untargeted channels remain."
            )

        for channel_name in unique_channels:
            current_sum = current_channel_sums[channel_name]
            desired_sum = class_total * desired_channel_shares[channel_name]
            if desired_sum > 0.0 and current_sum <= 0.0:
                raise RuntimeError(
                    f"Cannot assign non-zero target share to channel '{channel_name}' in class '{class_name}' because its current weight is zero."
                )
            scale = 0.0 if current_sum <= 0.0 else desired_sum / current_sum
            channel_mask = class_mask & (channel_names == channel_name)
            effective[channel_mask] *= np.float32(scale)
            summary.append(
                {
                    "class_name": class_name,
                    "channel_name": channel_name,
                    "before_sum": current_sum,
                    "after_sum": float(np.sum(effective[channel_mask], dtype=np.float64)),
                    "target_share": desired_channel_shares[channel_name],
                    "effective_scale": scale,
                }
            )

    return effective, summary


def summarize_channel_reweighting(weights_before, weights_after, labels, channel_ids, channel_info, class_names):
    summary = []
    unique_channel_ids = np.unique(channel_ids.astype(np.int64))
    for channel_id in unique_channel_ids:
        info = channel_info.get(int(channel_id), {})
        channel_mask = channel_ids == channel_id
        label_index = int(info.get("label_index", -1))
        label_name = class_names[label_index] if 0 <= label_index < len(class_names) else "unknown"
        before_sum = float(np.sum(weights_before[channel_mask], dtype=np.float64))
        after_sum = float(np.sum(weights_after[channel_mask], dtype=np.float64))
        summary.append(
            {
                "channel_id": int(channel_id),
                "original_name": info.get("original_name", f"channel_{int(channel_id)}"),
                "channel_name": info.get("channel_name", "unknown"),
                "class_name": info.get("class_name", "unknown"),
                "label_name": label_name,
                "n_events": int(np.sum(channel_mask)),
                "before_sum": before_sum,
                "after_sum": after_sum,
                "effective_scale": (after_sum / before_sum) if abs(before_sum) > 0.0 else 0.0,
            }
        )
    return summary


def print_channel_reweighting_summary(title, summary, applied_rules, preserve_class_sums):
    print(title)
    for item in applied_rules:
        print(item)
    for item in summary:
        if "channel_id" in item:
            print(
                f"  channel {item['channel_id']} ({item['original_name']} | {item['channel_name']} | {item['label_name']}): "
                f"n={item['n_events']}, before={item['before_sum']:.6g}, after={item['after_sum']:.6g}, "
                f"effective_scale={item['effective_scale']:.6g}"
            )
        else:
            print(
                f"  class {item['class_name']} / channel {item['channel_name']}: "
                f"before={item['before_sum']:.6g}, after={item['after_sum']:.6g}, "
                f"target_share={item['target_share']:.6g}, effective_scale={item['effective_scale']:.6g}"
            )


def print_class_target_summary(title, summary):
    print(title)
    for item in summary:
        print(
            f"  class {item['class_name']}: before={item['before_sum']:.6g}, after={item['after_sum']:.6g}, "
            f"target_relative={item['target_relative']:.6g}, effective_scale={item['effective_scale']:.6g}"
        )


def collect_weight_summary(weights, labels, class_names=None):
    if class_names is None:
        class_names = CLASS_NAMES
    summary = []
    for class_index, class_name in enumerate(class_names):
        class_mask = labels == class_index
        class_weights = weights[class_mask].astype(np.float64)
        event_count = int(np.sum(class_mask))
        weight_sum = float(np.sum(class_weights))
        mean_weight = float(np.mean(class_weights)) if event_count > 0 else 0.0
        min_weight = float(np.min(class_weights)) if event_count > 0 else 0.0
        max_weight = float(np.max(class_weights)) if event_count > 0 else 0.0
        summary.append(
            {
                "class_index": class_index,
                "class_name": class_name,
                "event_count": event_count,
                "weight_sum": weight_sum,
                "mean_weight": mean_weight,
                "min_weight": min_weight,
                "max_weight": max_weight,
            }
        )
    return summary


def print_weight_summary(title, weights, labels, class_names=None):
    print(title)
    for item in collect_weight_summary(weights, labels, class_names=class_names):
        print(
            f"  class {item['class_index']} ({item['class_name']}): "
            f"n={item['event_count']}, "
            f"sumW={item['weight_sum']:.6g}, "
            f"meanW={item['mean_weight']:.6g}, "
            f"minW={item['min_weight']:.6g}, "
            f"maxW={item['max_weight']:.6g}"
        )


def plot_weight_summary(raw_weights, normalized_weights, labels, output_path, title, class_names=None):
    if class_names is None:
        class_names = CLASS_NAMES
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    raw_summary = collect_weight_summary(raw_weights, labels, class_names=class_names)
    normalized_summary = collect_weight_summary(normalized_weights, labels, class_names=class_names)
    x = np.arange(len(class_names))
    width = 0.38

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))

    raw_sums = [item["weight_sum"] for item in raw_summary]
    norm_sums = [item["weight_sum"] for item in normalized_summary]
    axes[0].bar(x - width / 2, raw_sums, width=width, label="Raw")
    axes[0].bar(x + width / 2, norm_sums, width=width, label="Class-normalized")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(class_names, rotation=15, ha="right")
    axes[0].set_ylabel("Total class weight")
    axes[0].set_title("Per-class weight sums")
    axes[0].grid(alpha=0.25, axis="y")
    axes[0].legend()

    raw_means = [item["mean_weight"] for item in raw_summary]
    norm_means = [item["mean_weight"] for item in normalized_summary]
    axes[1].bar(x - width / 2, raw_means, width=width, label="Raw")
    axes[1].bar(x + width / 2, norm_means, width=width, label="Class-normalized")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(class_names, rotation=15, ha="right")
    axes[1].set_ylabel("Mean event weight")
    axes[1].set_title("Per-class mean weights")
    axes[1].grid(alpha=0.25, axis="y")
    axes[1].legend()

    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def weighted_cross_entropy(logits, targets, sample_weights):
    losses = nn.functional.cross_entropy(logits, targets, reduction="none")
    weight_sum = torch.clamp(sample_weights.sum(), min=1e-12)
    return torch.sum(losses * sample_weights) / weight_sum


def add_dataset_report_argument(parser):
    parser.add_argument(
        "--dataset-report",
        default=None,
        help=(
            f"Optional {DEFAULT_MIX_REPORT_NAME} from prepareSignalDatasets.py. "
            "If omitted, it is resolved next to the dataset file."
        ),
    )


def add_label_subset_argument(parser):
    parser.add_argument(
        "--label-subset",
        nargs="+",
        default=None,
        help=(
            "Optional subset of dataset labels used for a sub-task. "
            "Each entry can be a class index or class name from mix_report.json, e.g. '--label-subset 0 1' or '--label-subset TTHH tt_b'."
        ),
    )


def add_feature_group_argument(parser):
    parser.add_argument(
        "--feature-groups",
        nargs="+",
        default=None,
        help=(
            "Optional physics-group subset of inputs to keep. "
            f"Available groups: {', '.join(get_group_names())}. "
            "Numeric aliases 1..7 are also accepted."
        ),
    )


def resolve_train_outputs(args):
    if args.output_dir is None:
        raise RuntimeError("Provide --output-dir.")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    return {
        "model_path": output_dir / "model.pt",
        "loss_path": output_dir / "loss.png",
        "weight_summary_prefix": None if not args.save_weight_summary else output_dir / "weight_summary",
        "output_dir": output_dir,
    }


def resolve_test_outputs(args):
    prediction_path = Path(args.save_predictions) if args.save_predictions is not None else None
    plot_dir = Path(args.plot_dir) if args.plot_dir is not None else None
    confusion_prefix = None
    if args.output_dir is not None:
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        if prediction_path is None:
            prediction_path = output_dir / "predictions.npy"
        if plot_dir is None and args.save_basic_plots:
            plot_dir = output_dir / "score_plots"
        confusion_prefix = output_dir / "confusion_matrix"
    return {
        "prediction_path": prediction_path,
        "plot_dir": plot_dir,
        "confusion_prefix": confusion_prefix,
    }


def load_dataset_contract(dataset_path, dataset_report_path=None):
    report_path = resolve_dataset_report_path(dataset_path, explicit_report_path=dataset_report_path)
    report = load_json(report_path)
    validate_dataset_report(report, str(report_path))
    return report_path, report


def build_base_data_contract(dataset_report):
    return {
        "feature_names": list(dataset_report["feature_names"]),
        "mixed_field_names": list(dataset_report["mixed_field_names"]),
        "source_label_names": list(dataset_report["label_names"]),
    }


def build_task_data_contract(dataset_report, label_subset_indices, feature_groups=None):
    base = build_base_data_contract(dataset_report)
    source_label_names = base["source_label_names"]
    if label_subset_indices is None:
        label_subset_indices = list(range(len(source_label_names)))
    base["task_label_indices"] = [int(index) for index in label_subset_indices]
    base["task_label_names"] = [source_label_names[int(index)] for index in label_subset_indices]
    base["active_feature_groups"] = resolve_feature_groups(feature_groups)
    base["feature_layout"] = "compact"
    base["active_feature_names"] = get_active_feature_names(base["active_feature_groups"])
    return base


def normalize_data_contract(contract):
    if contract is None:
        return None

    normalized = dict(contract)
    if "source_label_names" not in normalized:
        label_names = list(normalized.get("label_names", []))
        normalized["source_label_names"] = label_names
    else:
        normalized["source_label_names"] = list(normalized["source_label_names"])

    if "task_label_indices" not in normalized or normalized["task_label_indices"] is None:
        normalized["task_label_indices"] = list(range(len(normalized["source_label_names"])))
    else:
        normalized["task_label_indices"] = [int(index) for index in normalized["task_label_indices"]]

    if "task_label_names" not in normalized or normalized["task_label_names"] is None:
        normalized["task_label_names"] = [
            normalized["source_label_names"][int(index)] for index in normalized["task_label_indices"]
        ]
    else:
        normalized["task_label_names"] = list(normalized["task_label_names"])

    normalized["feature_names"] = list(normalized["feature_names"])
    normalized["mixed_field_names"] = list(normalized["mixed_field_names"])
    if "active_feature_groups" not in normalized or normalized["active_feature_groups"] is None:
        normalized["active_feature_groups"] = list(get_group_names())
    else:
        normalized["active_feature_groups"] = resolve_feature_groups(normalized["active_feature_groups"])
    if "feature_layout" not in normalized or normalized["feature_layout"] is None:
        normalized["feature_layout"] = "masked"
    if "active_feature_names" not in normalized or normalized["active_feature_names"] is None:
        normalized["active_feature_names"] = get_active_feature_names(normalized["active_feature_groups"])
    else:
        normalized["active_feature_names"] = list(normalized["active_feature_names"])
    return normalized


def resolve_label_subset(label_subset_args, source_label_names):
    if label_subset_args is None:
        return None

    resolved = []
    for token in label_subset_args:
        if token.isdigit():
            index = int(token)
            if index < 0 or index >= len(source_label_names):
                raise RuntimeError(f"Label index {index} is outside [0, {len(source_label_names) - 1}].")
        else:
            if token not in source_label_names:
                raise RuntimeError(
                    f"Unknown label name '{token}'. Available: {', '.join(source_label_names)}"
                )
            index = source_label_names.index(token)
        if index not in resolved:
            resolved.append(index)

    if len(resolved) < 2:
        raise RuntimeError("label-subset must contain at least two distinct labels.")
    return resolved


def filter_dataset_for_task(x, weights, labels, is_same_sign, channel_ids, label_subset_indices):
    if label_subset_indices is None:
        return x, weights, labels, is_same_sign, channel_ids

    keep_mask = build_task_keep_mask(labels, label_subset_indices)
    if not np.any(keep_mask):
        raise RuntimeError("No events remain after applying --label-subset.")

    subset = np.asarray(label_subset_indices, dtype=np.int64)
    remap = {int(label): new_index for new_index, label in enumerate(subset.tolist())}
    filtered_labels = labels[keep_mask].astype(np.int64).copy()
    remapped_labels = np.array([remap[int(label)] for label in filtered_labels], dtype=np.int64)
    return (
        x[keep_mask],
        weights[keep_mask],
        remapped_labels,
        is_same_sign[keep_mask],
        channel_ids[keep_mask],
    )


def build_task_keep_mask(labels, label_subset_indices):
    if label_subset_indices is None:
        return np.ones(labels.shape[0], dtype=bool)
    subset = np.asarray(label_subset_indices, dtype=np.int64)
    return np.isin(labels, subset)


def save_model_checkpoint(model, output_path, data_contract):
    checkpoint = {
        "state_dict": model.state_dict(),
        "data_contract": data_contract,
    }
    torch.save(checkpoint, output_path)


def validate_checkpoint_contract(saved_contract, requested_contract, model_path):
    if requested_contract is None or saved_contract is None:
        return
    compare_task = "task_label_indices" in requested_contract
    saved_contract = normalize_data_contract(saved_contract)
    requested_contract = normalize_data_contract(requested_contract)
    keys_to_match = ("feature_names", "mixed_field_names", "source_label_names")
    for key in keys_to_match:
        if saved_contract.get(key) != requested_contract.get(key):
            raise RuntimeError(f"{model_path}: checkpoint {key} does not match the dataset report.")
    if compare_task:
        if saved_contract.get("task_label_indices") != requested_contract.get("task_label_indices"):
            raise RuntimeError(f"{model_path}: checkpoint task_label_indices do not match the requested sub-task.")
        if saved_contract.get("task_label_names") != requested_contract.get("task_label_names"):
            raise RuntimeError(f"{model_path}: checkpoint task_label_names do not match the requested sub-task.")
    if saved_contract.get("active_feature_groups") != requested_contract.get("active_feature_groups"):
        raise RuntimeError(f"{model_path}: checkpoint active_feature_groups do not match the requested feature mask.")
    if saved_contract.get("feature_layout") != requested_contract.get("feature_layout"):
        raise RuntimeError(f"{model_path}: checkpoint feature_layout does not match the requested input layout.")
    if saved_contract.get("active_feature_names") != requested_contract.get("active_feature_names"):
        raise RuntimeError(f"{model_path}: checkpoint active_feature_names do not match the requested features.")


def load_model_checkpoint(model_path, device, requested_contract=None, requested_schema=None):
    if requested_contract is None and requested_schema is not None:
        requested_contract = {
            "feature_names": list(FEATURE_FIELD_NAMES),
            "mixed_field_names": list(MIXED_FIELD_NAMES),
            "source_label_names": list(get_signal_schema(requested_schema).class_names),
        }
    checkpoint = torch.load(model_path, map_location=device)
    if isinstance(checkpoint, dict) and "state_dict" in checkpoint:
        state_dict = checkpoint["state_dict"]
        saved_contract = checkpoint.get("data_contract")
        if saved_contract is None:
            saved_schema = checkpoint.get("class_schema", DEFAULT_SCHEMA_NAME)
            saved_contract = {
                "feature_names": list(FEATURE_FIELD_NAMES),
                "mixed_field_names": list(MIXED_FIELD_NAMES),
                "source_label_names": list(get_signal_schema(saved_schema).class_names),
            }
    else:
        state_dict = checkpoint
        saved_contract = requested_contract

    saved_contract = normalize_data_contract(saved_contract)
    validate_checkpoint_contract(saved_contract, requested_contract, model_path)
    return state_dict, saved_contract


def parse_args():
    parser = argparse.ArgumentParser(description="Train or evaluate a simple signal classification MLP.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    train_parser = subparsers.add_parser("train")
    train_parser.add_argument("--train", required=True, help="Training .npy dataset.")
    train_parser.add_argument("--val", required=True, help="Validation .npy dataset.")
    train_parser.add_argument("--output-dir", required=True, help="Directory for training outputs.")
    train_parser.add_argument("--epochs", type=int, default=100)
    train_parser.add_argument("--batch-size", type=int, default=1024)
    train_parser.add_argument("--lr", type=float, default=1e-3)
    train_parser.add_argument("--device", default=None)
    train_parser.add_argument(
        "--save-weight-summary",
        action="store_true",
        help="Optional: also save train/val weight summary plots.",
    )
    train_parser.add_argument(
        "--training-weight-mode",
        choices=["balanced", "input", "abs"],
        default="balanced",
        help=(
            "How to build training weights from the dataset weight column. "
            "'balanced' keeps the previous behavior; 'input' uses the dataset weights as-is; "
            "'abs' uses abs(weight) without class re-normalization."
        ),
    )
    train_parser.add_argument(
        "--weight-target",
        action="append",
        default=[],
        help=(
            "Unified weight interface. "
            "Use 'class=value' for inter-class relative importance, e.g. '--weight-target tt_b=2'. "
            "Use 'class/channel=share' for in-class fractions, e.g. '--weight-target tt_b/TTBB=0.8'. "
            "Repeat as needed."
        ),
    )
    add_dataset_report_argument(train_parser)
    add_label_subset_argument(train_parser)
    add_feature_group_argument(train_parser)
    train_parser.add_argument(
        "--feature-layout",
        choices=["compact", "masked"],
        default="compact",
        help=(
            "'compact' keeps only active feature columns and reduces model input dimension. "
            "'masked' keeps the original full input dimension and zeros inactive columns."
        ),
    )

    test_parser = subparsers.add_parser("test")
    test_parser.add_argument("--test", required=True, help="Test .npy dataset.")
    test_parser.add_argument("--model", required=True, help="Trained model .pt file.")
    test_parser.add_argument("--output-dir", default=None, help="Directory for test outputs.")
    test_parser.add_argument("--batch-size", type=int, default=1024)
    test_parser.add_argument("--device", default=None)
    test_parser.add_argument("--suffix", default="_simple_signal")
    test_parser.add_argument("--save-predictions", default=None, help="Optional output .npy for structured predictions.")
    test_parser.add_argument("--plot-dir", default=None, help="Optional output directory for score histograms.")
    test_parser.add_argument(
        "--save-basic-plots",
        action="store_true",
        help="If set, save the basic score histograms into --output-dir/score_plots when --plot-dir is not given.",
    )
    test_parser.add_argument(
        "--score-bins",
        nargs="+",
        type=float,
        default=[0.0, 0.2, 0.5, 0.8, 1.0],
        help="Histogram bin edges for score plots.",
    )
    test_parser.add_argument(
        "--score-normalization",
        choices=["raw", "global", "shape", "shape_abs"],
        default="shape",
        help=(
            "Score-plot normalization mode. "
            "'raw' keeps original event weights, "
            "'global' divides all bins by the total weight across the plotted sample, "
            "'shape' normalizes each truth class separately using signed weights, "
            "'shape_abs' normalizes each truth class separately using abs(weight)."
        ),
    )
    test_parser.add_argument(
        "--equalize-classes",
        action="store_true",
        help="Deprecated alias for '--score-normalization shape'.",
    )
    add_dataset_report_argument(test_parser)
    add_label_subset_argument(test_parser)
    add_feature_group_argument(test_parser)
    test_parser.add_argument(
        "--feature-layout",
        choices=["compact", "masked"],
        default=None,
        help="Optional override for the feature layout. Defaults to the layout stored in the checkpoint.",
    )

    return parser.parse_args()


def apply_feature_layout(features, active_feature_groups, layout):
    if layout == "masked":
        return apply_feature_mask(features, active_feature_groups)
    if layout == "compact":
        return select_active_features(features, active_feature_groups)
    raise RuntimeError(f"Unknown feature layout: {layout}")


def get_cuda_diagnostics():
    diagnostics = [
        f"torch={torch.__version__}",
        f"torch_cuda={torch.version.cuda}",
    ]
    try:
        with warnings.catch_warnings(record=True) as caught_warnings:
            warnings.simplefilter("always")
            is_available = torch.cuda.is_available()
        diagnostics.append(f"cuda_available={is_available}")
        if is_available:
            diagnostics.append(f"cuda_device_count={torch.cuda.device_count()}")
            diagnostics.append(f"cuda_device_name={torch.cuda.get_device_name(0)}")
        elif caught_warnings:
            diagnostics.append(f"cuda_warning={caught_warnings[-1].message}")
    except Exception as exc:
        diagnostics.append(f"cuda_error={exc}")
    return " | ".join(str(item) for item in diagnostics)


def resolve_device(device_name):
    if device_name is not None:
        device = torch.device(device_name)
        if device.type == "cuda":
            if not torch.cuda.is_available():
                raise RuntimeError(
                    f"Requested CUDA device '{device_name}', but CUDA is not available. {get_cuda_diagnostics()}"
                )
            print(f"Using device: {device} | {get_cuda_diagnostics()}")
        else:
            print(f"Using device: {device}")
        return device

    if torch.cuda.is_available():
        device = torch.device("cuda")
        print(f"Using device: {device} | {get_cuda_diagnostics()}")
        return device

    print(f"Using device: cpu | {get_cuda_diagnostics()}")
    return torch.device("cpu")


def train_model(args):
    report_path, dataset_report = load_dataset_contract(args.train, args.dataset_report)
    channel_info = load_channel_info_from_report(dataset_report)
    source_label_names = list(dataset_report["label_names"])
    label_subset_indices = resolve_label_subset(args.label_subset, source_label_names)
    feature_groups = resolve_feature_groups(args.feature_groups)
    task_contract = build_task_data_contract(dataset_report, label_subset_indices, feature_groups=feature_groups)
    task_contract["feature_layout"] = args.feature_layout
    task_contract["training_weight_mode"] = args.training_weight_mode
    class_names = list(task_contract["task_label_names"])
    validate_mixed_dataset_file(args.train, dataset_report, expected_sign="os")
    validate_mixed_dataset_file(args.val, dataset_report, expected_sign="os")
    print(f"validated dataset report: {report_path}")
    outputs = resolve_train_outputs(args)
    x_train, train_weights_raw, y_train, _, train_channel_ids = loadDS(args.train)
    x_val, val_weights_raw, y_val, _, val_channel_ids = loadDS(args.val)
    x_train, train_weights_raw, y_train, _, train_channel_ids = filter_dataset_for_task(
        x_train, train_weights_raw, y_train, np.zeros_like(y_train), train_channel_ids, label_subset_indices
    )
    x_val, val_weights_raw, y_val, _, val_channel_ids = filter_dataset_for_task(
        x_val, val_weights_raw, y_val, np.zeros_like(y_val), val_channel_ids, label_subset_indices
    )
    x_train = apply_feature_layout(x_train, feature_groups, args.feature_layout)
    x_val = apply_feature_layout(x_val, feature_groups, args.feature_layout)
    train_weights_abs = make_non_negative_training_weights(train_weights_raw)
    val_weights_abs = make_non_negative_training_weights(val_weights_raw)
    train_weights = resolve_training_weights(train_weights_raw, y_train, args.training_weight_mode)
    val_weights = resolve_training_weights(val_weights_raw, y_val, args.training_weight_mode)
    class_targets, channel_share_targets = resolve_weight_targets(args)
    train_weights_before_class = train_weights.copy()
    val_weights_before_class = val_weights.copy()
    train_weights, train_class_target_summary = apply_class_weight_targets(
        train_weights,
        y_train,
        class_names,
        class_targets,
    )
    val_weights, val_class_target_summary = apply_class_weight_targets(
        val_weights,
        y_val,
        class_names,
        class_targets,
    )
    train_weights_before_channel = train_weights.copy()
    val_weights_before_channel = val_weights.copy()
    train_weights, train_channel_target_summary = apply_class_channel_share_targets(
        train_weights,
        y_train,
        train_channel_ids,
        channel_info,
        class_names,
        channel_share_targets,
    )
    val_weights, val_channel_target_summary = apply_class_channel_share_targets(
        val_weights,
        y_val,
        val_channel_ids,
        channel_info,
        class_names,
        channel_share_targets,
    )
    task_contract["weight_targets"] = {
        "class_targets": class_targets,
        "channel_share_targets": channel_share_targets,
    }

    print_weight_summary("Train raw weights:", train_weights_raw, y_train, class_names=class_names)
    print_weight_summary("Train abs(weights):", train_weights_abs, y_train, class_names=class_names)
    print_weight_summary(
        f"Train base weights before explicit targets ({args.training_weight_mode}):",
        train_weights_before_class,
        y_train,
        class_names=class_names,
    )
    print_class_target_summary("Train class-target summary:", train_class_target_summary)
    print_channel_reweighting_summary(
        "Train class/channel-target summary:",
        train_channel_target_summary,
        [],
        None,
    )
    print_weight_summary(
        f"Train effective training weights ({args.training_weight_mode}):",
        train_weights,
        y_train,
        class_names=class_names,
    )
    print_weight_summary("Val raw weights:", val_weights_raw, y_val, class_names=class_names)
    print_weight_summary("Val abs(weights):", val_weights_abs, y_val, class_names=class_names)
    print_weight_summary(
        f"Val base weights before explicit targets ({args.training_weight_mode}):",
        val_weights_before_class,
        y_val,
        class_names=class_names,
    )
    print_class_target_summary("Val class-target summary:", val_class_target_summary)
    print_channel_reweighting_summary(
        "Val class/channel-target summary:",
        val_channel_target_summary,
        [],
        None,
    )
    print_weight_summary(
        f"Val effective training weights ({args.training_weight_mode}):",
        val_weights,
        y_val,
        class_names=class_names,
    )
    print_channel_reweighting_summary(
        "Train final per-channel effect:",
        summarize_channel_reweighting(
            train_weights_before_channel,
            train_weights,
            y_train,
            train_channel_ids,
            channel_info,
            class_names,
        ),
        [],
        None,
    )
    print_channel_reweighting_summary(
        "Val final per-channel effect:",
        summarize_channel_reweighting(
            val_weights_before_channel,
            val_weights,
            y_val,
            val_channel_ids,
            channel_info,
            class_names,
        ),
        [],
        None,
    )

    if outputs["weight_summary_prefix"] is not None:
        weight_summary_prefix = outputs["weight_summary_prefix"]
        plot_weight_summary(
            train_weights_raw,
            train_weights,
            y_train,
            Path(str(weight_summary_prefix) + "_train.png"),
            f"Train weight summary ({args.training_weight_mode})",
            class_names=class_names,
        )
        plot_weight_summary(
            val_weights_raw,
            val_weights,
            y_val,
            Path(str(weight_summary_prefix) + "_val.png"),
            f"Validation weight summary ({args.training_weight_mode})",
            class_names=class_names,
        )

    model = SimpleMLP(input_dim=x_train.shape[1], n_classes=len(class_names))
    print_training_start_summary(args, model, class_names, report_path)
    device = resolve_device(args.device)
    model = model.to(device)

    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)

    train_loader = build_loader(x_train, y_train, train_weights, args.batch_size, shuffle=True)
    val_loader = build_loader(x_val, y_val, val_weights, args.batch_size, shuffle=False)

    train_losses = []
    val_losses = []
    for epoch in range(args.epochs):
        model.train()
        total_train_loss = 0.0
        for batch_x, batch_y, batch_w in train_loader:
            batch_x = batch_x.to(device)
            batch_y = batch_y.to(device)
            batch_w = batch_w.to(device)

            optimizer.zero_grad()
            logits = model(batch_x)
            loss = weighted_cross_entropy(logits, batch_y, batch_w)
            loss.backward()
            optimizer.step()
            total_train_loss += loss.item()

        train_losses.append(total_train_loss / max(len(train_loader), 1))

        model.eval()
        total_val_loss = 0.0
        correct = 0
        with torch.no_grad():
            for batch_x, batch_y, batch_w in val_loader:
                batch_x = batch_x.to(device)
                batch_y = batch_y.to(device)
                batch_w = batch_w.to(device)
                logits = model(batch_x)
                loss = weighted_cross_entropy(logits, batch_y, batch_w)
                total_val_loss += loss.item()
                correct += (logits.argmax(dim=1) == batch_y).sum().item()

        val_losses.append(total_val_loss / max(len(val_loader), 1))
        val_acc = correct / max(len(y_val), 1)
        print(
            f"Epoch {epoch + 1}/{args.epochs} | "
            f"Train Loss: {train_losses[-1]:.4f} | "
            f"Val Loss: {val_losses[-1]:.4f} | "
            f"Val Acc: {val_acc:.4f}"
        )

    save_model_checkpoint(model, outputs["model_path"], task_contract)
    plot_loss_curves(train_losses, val_losses, outputs["loss_path"])
    print(f"saved model: {outputs['model_path']}")
    print(f"saved loss plot: {outputs['loss_path']}")
    print(f"training weight mode: {args.training_weight_mode}")


def test_model(args):
    outputs = resolve_test_outputs(args)
    report_path, dataset_report = load_dataset_contract(args.test, args.dataset_report)
    validate_mixed_dataset_file(args.test, dataset_report)
    base_contract = build_base_data_contract(dataset_report)
    source_label_names = base_contract["source_label_names"]
    requested_subset = resolve_label_subset(args.label_subset, source_label_names)
    requested_feature_groups = resolve_feature_groups(args.feature_groups)
    device = resolve_device(args.device)
    state_dict, checkpoint_contract = load_model_checkpoint(args.model, device)
    effective_subset = checkpoint_contract["task_label_indices"] if requested_subset is None else requested_subset
    effective_feature_groups = (
        checkpoint_contract["active_feature_groups"] if args.feature_groups is None else requested_feature_groups
    )
    effective_feature_layout = checkpoint_contract["feature_layout"] if args.feature_layout is None else args.feature_layout
    requested_contract = build_task_data_contract(
        dataset_report,
        effective_subset,
        feature_groups=effective_feature_groups,
    )
    requested_contract["feature_layout"] = effective_feature_layout
    validate_checkpoint_contract(checkpoint_contract, requested_contract, args.model)
    class_names = list(requested_contract["task_label_names"])
    x_test, _, y_test, is_same_sign, channel_ids = loadDS(args.test)
    task_keep_mask = build_task_keep_mask(y_test, effective_subset)
    x_test, _, y_test, is_same_sign, channel_ids = filter_dataset_for_task(
        x_test, np.ones_like(y_test, dtype=np.float32), y_test, is_same_sign, channel_ids, effective_subset
    )
    x_test = apply_feature_layout(x_test, effective_feature_groups, effective_feature_layout)
    model = SimpleMLP(input_dim=x_test.shape[1], n_classes=len(class_names))
    model.load_state_dict(state_dict)
    model = model.to(device)
    model.eval()

    test_loader = build_loader(x_test, y_test, None, args.batch_size, shuffle=False)
    prediction = []
    truth = []
    all_scores = []
    with torch.no_grad():
        for batch_x, batch_y in test_loader:
            batch_x = batch_x.to(device)
            logits = model(batch_x)
            score_batch = torch.softmax(logits, dim=1).cpu().numpy()
            all_scores.append(score_batch)
            prediction.extend(np.argmax(score_batch, axis=1).tolist())
            truth.extend(batch_y.numpy().tolist())

    prediction = np.array(prediction, dtype=np.int64)
    truth = np.array(truth, dtype=np.int64)
    scores = np.concatenate(all_scores, axis=0).astype(np.float32) if all_scores else np.empty((0, len(class_names)), dtype=np.float32)
    accuracy = np.mean(prediction == truth)
    print("acc:", accuracy)
    print(f"validated dataset report: {report_path}")
    print("task labels:", ", ".join(class_names))
    if is_same_sign.size > 0:
        print("same-sign fraction:", np.mean(is_same_sign))
    plot_confusion_matrix_fractional(
        truth,
        prediction,
        args.suffix,
        class_names,
        class_names,
        output_name_prefix="confusion_matrix" if outputs["confusion_prefix"] is None else str(outputs["confusion_prefix"]),
    )
    if outputs["prediction_path"] is not None:
        raw_dataset = np.load(args.test)
        raw_matrix = structured_to_matrix(raw_dataset)
        raw_weights = raw_matrix[:, WEIGHT_INDEX].astype(np.float32)
        _, filtered_weights, _, filtered_signs, _ = filter_dataset_for_task(
            raw_matrix[:, :FEATURE_DIM].astype(np.float32),
            raw_weights,
            raw_matrix[:, LABEL_INDEX].astype(np.int64),
            raw_matrix[:, SIGN_INDEX].astype(np.float32),
            raw_matrix[:, CHANNEL_INDEX].astype(np.int64),
            effective_subset,
        )
        original_row_index = np.flatnonzero(task_keep_mask).astype(np.int64)
        save_structured_predictions(
            outputs["prediction_path"],
            truth,
            prediction,
            scores,
            filtered_weights,
            filtered_signs,
            class_names,
            row_index=original_row_index,
        )
        print(f"saved predictions: {outputs['prediction_path']}")
    if outputs["plot_dir"] is not None:
        raw_dataset = np.load(args.test)
        raw_matrix = structured_to_matrix(raw_dataset)
        filtered_features, filtered_weights, _, _, _ = filter_dataset_for_task(
            raw_matrix[:, :FEATURE_DIM].astype(np.float32),
            raw_matrix[:, WEIGHT_INDEX].astype(np.float32),
            raw_matrix[:, LABEL_INDEX].astype(np.int64),
            raw_matrix[:, SIGN_INDEX].astype(np.float32),
            raw_matrix[:, CHANNEL_INDEX].astype(np.int64),
            effective_subset,
        )
        filtered_features = apply_feature_layout(filtered_features, effective_feature_groups, effective_feature_layout)
        score_normalization = "shape" if args.equalize_classes else args.score_normalization
        plot_score_histograms(
            truth=truth,
            scores=scores,
            weights=filtered_weights,
            output_dir=outputs["plot_dir"],
            suffix=args.suffix,
            bin_edges=np.asarray(args.score_bins, dtype=np.float32),
            normalization_mode=score_normalization,
            class_names=class_names,
        )
        print(f"saved basic score plots: {outputs['plot_dir']}")


def extract_event_weights(raw_dataset):
    matrix = structured_to_matrix(raw_dataset)
    return matrix[:, WEIGHT_INDEX].astype(np.float32)


def save_structured_predictions(output_path, truth, prediction, scores, weights, signs, class_names, row_index=None):
    if row_index is None:
        row_index = np.arange(len(truth), dtype=np.int64)
    dtype = [
        ("rowIndex", np.int64),
        ("truth", np.int32),
        ("prediction", np.int32),
        ("weight", np.float32),
        ("isSameSign", np.float32),
    ]
    dtype.extend((class_name, np.float32) for class_name in class_names)
    output = np.empty(len(truth), dtype=dtype)
    output["rowIndex"] = row_index
    output["truth"] = truth
    output["prediction"] = prediction
    output["weight"] = weights
    output["isSameSign"] = signs
    for class_index, class_name in enumerate(class_names):
        output[class_name] = scores[:, class_index]
    np.save(output_path, output)


def plot_score_histograms(truth, scores, weights, output_dir, suffix, bin_edges, normalization_mode, class_names):
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    palette = plt.get_cmap("tab10").colors
    total_weight = np.sum(weights.astype(np.float64), dtype=np.float64)

    for score_index, score_name in enumerate(class_names):
        fig, ax = plt.subplots(figsize=(8, 6))
        for class_index, class_label in enumerate(class_names):
            class_mask = truth == class_index
            class_scores = scores[class_mask, score_index]
            class_weights = weights[class_mask].astype(np.float64)
            if normalization_mode == "global":
                if np.abs(total_weight) > 0:
                    class_weights = class_weights / total_weight
            elif normalization_mode == "shape":
                class_weight_sum = np.sum(class_weights, dtype=np.float64)
                if np.abs(class_weight_sum) > 0:
                    class_weights = class_weights / class_weight_sum
            elif normalization_mode == "shape_abs":
                class_weights = np.abs(class_weights)
                class_weight_sum = np.sum(class_weights, dtype=np.float64)
                if class_weight_sum > 0:
                    class_weights = class_weights / class_weight_sum
            ax.hist(
                class_scores,
                bins=bin_edges,
                weights=class_weights,
                histtype="step",
                linewidth=2.0,
                color=palette[class_index % len(palette)],
                label=class_label,
            )

        ax.set_xlabel(f"{score_name} score")
        if normalization_mode == "raw":
            ax.set_ylabel("Weighted entries")
        elif normalization_mode == "global":
            ax.set_ylabel("Globally normalized weighted entries")
        elif normalization_mode == "shape":
            ax.set_ylabel("Per-class normalized weighted entries")
        else:
            ax.set_ylabel("Per-class normalized abs-weighted entries")
        ax.set_title(f"{score_name} score by truth class")
        ax.set_xlim(float(bin_edges[0]), float(bin_edges[-1]))
        ax.set_xticks(bin_edges)
        ax.axhline(0.0, color="black", linewidth=0.8, alpha=0.5)
        ax.grid(alpha=0.25)
        ax.legend(title="Truth class")
        fig.tight_layout()
        fig.savefig(output_path / f"{score_name}{suffix}.png")
        plt.close(fig)


def main():
    args = parse_args()
    if args.command == "train":
        train_model(args)
    elif args.command == "test":
        test_model(args)


if __name__ == "__main__":
    main()
