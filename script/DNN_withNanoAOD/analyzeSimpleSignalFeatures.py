import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import torch
from torch.utils.data import DataLoader, TensorDataset

from feature_grouping import apply_feature_mask, build_feature_groups, select_active_features
from signalExtractionInputFromNanoAOD import FEATURE_FIELD_NAMES
from prepareSignalDatasets import (
    DEFAULT_MIX_REPORT_NAME,
    load_json,
    resolve_dataset_report_path,
    validate_dataset_report,
    validate_mixed_dataset_file,
)
from simpleSignalModel import (
    SimpleMLP,
    build_balanced_training_weights,
    build_base_data_contract,
    build_task_data_contract,
    filter_dataset_for_task,
    loadDS,
    load_model_checkpoint,
    resolve_device,
    resolve_label_subset,
    validate_checkpoint_contract,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Analyze feature behavior for the simple signal model."
    )
    parser.add_argument("--dataset", required=True, help="Input dataset .npy, e.g. OS_train.npy.")
    parser.add_argument("--model", required=True, help="Trained .pt model.")
    parser.add_argument("--output-dir", default="feature_analysis", help="Output directory.")
    parser.add_argument(
        "--dataset-report",
        default=None,
        help=(
            f"Optional {DEFAULT_MIX_REPORT_NAME} from prepareSignalDatasets.py. "
            "If omitted, it is resolved next to the dataset file."
        ),
    )
    parser.add_argument("--batch-size", type=int, default=1024, help="Batch size for gradient analysis.")
    parser.add_argument("--device", default=None, help="cpu or cuda")
    parser.add_argument(
        "--label-subset",
        nargs="+",
        default=None,
        help="Optional sub-task labels, e.g. '--label-subset 0 1' or '--label-subset TTHH tt_b'.",
    )
    parser.add_argument(
        "--weight-mode",
        choices=["balanced", "abs", "raw", "unit"],
        default="balanced",
        help="Sample weighting used in gradient/covariance summaries.",
    )
    parser.add_argument("--perm-repeats", type=int, default=5, help="Permutation repeats per feature.")
    parser.add_argument("--perm-max-events", type=int, default=20000, help="Max events used for permutation importance.")
    parser.add_argument("--perm-min-events", type=int, default=200, help="Minimum events required to run a channel-wise permutation study.")
    parser.add_argument("--seed", type=int, default=1234, help="Random seed.")
    parser.add_argument(
        "--gradient-views",
        choices=["basic", "full"],
        default="basic",
        help="'basic' saves only the by-target gradient view. 'full' also saves overall and true-class gradient views.",
    )
    return parser.parse_args()


def resolve_analysis_weights(raw_weights, labels, mode):
    if mode == "balanced":
        return build_balanced_training_weights(raw_weights, labels).astype(np.float64)
    if mode == "abs":
        return np.abs(raw_weights).astype(np.float64)
    if mode == "raw":
        return raw_weights.astype(np.float64)
    if mode == "unit":
        return np.ones_like(raw_weights, dtype=np.float64)
    raise RuntimeError(f"Unknown weight mode: {mode}")


def load_model(model_path, device, class_names, data_contract):
    model = SimpleMLP(n_classes=len(class_names)).to(device)
    state_dict, _ = load_model_checkpoint(model_path, device, requested_contract=data_contract)
    model.load_state_dict(state_dict)
    model.eval()
    return model


def load_dataset_contract(dataset_path, dataset_report_path=None):
    report_path = resolve_dataset_report_path(dataset_path, explicit_report_path=dataset_report_path)
    report = load_json(report_path)
    validate_dataset_report(report, str(report_path))
    return report_path, report


def load_channel_id_to_info(dataset_report):
    raw_info = dataset_report["labeling"]["channel_info"]
    return {int(channel_id): info for channel_id, info in raw_info.items()}


def weighted_mean(values, weights):
    weight_sum = np.sum(weights, dtype=np.float64)
    if np.abs(weight_sum) == 0:
        return np.zeros(values.shape[1], dtype=np.float64)
    return np.sum(values * weights[:, None], axis=0, dtype=np.float64) / weight_sum


def weighted_covariance(values, weights):
    mean = weighted_mean(values, weights)
    centered = values - mean
    weight_sum = np.sum(weights, dtype=np.float64)
    if np.abs(weight_sum) == 0:
        return np.zeros((values.shape[1], values.shape[1]), dtype=np.float64)
    cov = (centered * weights[:, None]).T @ centered / weight_sum
    return cov.astype(np.float64)


def covariance_to_correlation(cov):
    diag = np.sqrt(np.maximum(np.diag(cov), 0.0))
    denom = np.outer(diag, diag)
    corr = np.zeros_like(cov)
    mask = denom > 0
    corr[mask] = cov[mask] / denom[mask]
    return np.clip(corr, -1.0, 1.0)


def compute_grad_importance(model, features, weights, batch_size, device, target_index=None):
    dataset = torch.utils.data.TensorDataset(
        torch.from_numpy(features).float(),
        torch.from_numpy(weights).float(),
    )
    loader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=False, drop_last=False)

    importance_sum = torch.zeros(features.shape[1], dtype=torch.float64, device=device)
    total_weight = torch.tensor(0.0, dtype=torch.float64, device=device)

    for batch_x, batch_w in loader:
        batch_x = batch_x.to(device)
        batch_w = batch_w.to(device)
        batch_x.requires_grad_(True)
        logits = model(batch_x)
        if target_index is None:
            target = logits.argmax(dim=1)
        else:
            target = torch.full((batch_x.shape[0],), int(target_index), dtype=torch.long, device=device)
        selected = logits[torch.arange(batch_x.shape[0], device=device), target].sum()
        grads = torch.autograd.grad(selected, batch_x, retain_graph=False, create_graph=False)[0]
        contribution = grads.abs() * batch_w[:, None]
        importance_sum += contribution.sum(dim=0).to(torch.float64)
        total_weight += batch_w.to(torch.float64).sum()

    if float(total_weight.cpu()) > 0.0:
        importance = importance_sum / total_weight
    else:
        importance = importance_sum
    return importance.detach().cpu().numpy().astype(np.float64)


def compute_binary_margin_importance(model, features, weights, batch_size, device):
    if features.shape[0] == 0:
        return np.zeros(features.shape[1], dtype=np.float64)

    dataset = torch.utils.data.TensorDataset(
        torch.from_numpy(features).float(),
        torch.from_numpy(weights).float(),
    )
    loader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=False, drop_last=False)

    importance_sum = torch.zeros(features.shape[1], dtype=torch.float64, device=device)
    total_weight = torch.tensor(0.0, dtype=torch.float64, device=device)

    for batch_x, batch_w in loader:
        batch_x = batch_x.to(device)
        batch_w = batch_w.to(device)
        batch_x.requires_grad_(True)
        logits = model(batch_x)
        if logits.shape[1] != 2:
            raise RuntimeError("compute_binary_margin_importance expects a binary classifier.")
        margin = (logits[:, 1] - logits[:, 0]).sum()
        grads = torch.autograd.grad(margin, batch_x, retain_graph=False, create_graph=False)[0]
        contribution = grads.abs() * batch_w[:, None]
        importance_sum += contribution.sum(dim=0).to(torch.float64)
        total_weight += batch_w.to(torch.float64).sum()

    if float(total_weight.cpu()) > 0.0:
        importance = importance_sum / total_weight
    else:
        importance = importance_sum
    return importance.detach().cpu().numpy().astype(np.float64)


@torch.no_grad()
def weighted_accuracy(model, features, labels, weights, batch_size, device):
    dataset = TensorDataset(
        torch.from_numpy(features).float(),
        torch.from_numpy(labels).long(),
        torch.from_numpy(weights).float(),
    )
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=False, drop_last=False)
    total_correct = 0.0
    total_weight = 0.0
    for batch_x, batch_y, batch_w in loader:
        batch_x = batch_x.to(device)
        batch_y = batch_y.to(device)
        logits = model(batch_x)
        pred = logits.argmax(dim=1)
        correct = (pred == batch_y).float() * batch_w.to(device)
        total_correct += float(correct.sum().cpu())
        total_weight += float(batch_w.sum())
    return total_correct / total_weight if total_weight > 0 else 0.0


@torch.no_grad()
def compute_permutation_importance(
    model, features, labels, weights, batch_size, device, n_repeats, seed, feature_indices=None
):
    rng = np.random.default_rng(seed)
    base_score = weighted_accuracy(model, features, labels, weights, batch_size, device)
    if feature_indices is None:
        feature_indices = list(range(features.shape[1]))
    importance = np.zeros(len(feature_indices), dtype=np.float64)
    work = features.copy()
    for output_index, feature_index in enumerate(feature_indices):
        drops = []
        original = work[:, feature_index].copy()
        for _ in range(n_repeats):
            perm = rng.permutation(work.shape[0])
            work[:, feature_index] = original[perm]
            score = weighted_accuracy(model, work, labels, weights, batch_size, device)
            drops.append(base_score - score)
        work[:, feature_index] = original
        importance[output_index] = max(0.0, float(np.mean(drops)))
    return importance, base_score


def compute_group_permutation_importance(
    model,
    features,
    labels,
    weights,
    batch_size,
    device,
    groups,
    n_repeats,
    seed,
):
    rng = np.random.default_rng(seed)
    base_score = weighted_accuracy(model, features, labels, weights, batch_size, device)
    importance = np.zeros(len(groups), dtype=np.float64)
    work = features.copy()
    for group_index, (_, feature_indices) in enumerate(groups):
        drops = []
        original = work[:, feature_indices].copy()
        for _ in range(n_repeats):
            perm = rng.permutation(work.shape[0])
            work[:, feature_indices] = original[perm]
            score = weighted_accuracy(model, work, labels, weights, batch_size, device)
            drops.append(base_score - score)
        work[:, feature_indices] = original
        importance[group_index] = max(0.0, float(np.mean(drops)))
    return importance, base_score


def plot_named_bar(values_by_label, names, output_path, title, ylabel):
    fig_height = max(6.0, 0.22 * len(names))
    fig, ax = plt.subplots(figsize=(11, fig_height))
    y_pos = np.arange(len(names))
    n_series = len(values_by_label)
    bar_height = 0.8 / max(n_series, 1)

    for idx, (label, values) in enumerate(values_by_label):
        offset = (idx - (n_series - 1) / 2.0) * bar_height
        ax.barh(y_pos + offset, values, height=bar_height, label=label, alpha=0.85)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(names, fontsize=8)
    ax.set_xlabel(ylabel)
    ax.set_title(title)
    ax.grid(alpha=0.25, axis="x")
    if n_series > 1:
        ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def plot_feature_bar(values_by_label, feature_names, output_path, title, ylabel):
    plot_named_bar(values_by_label, feature_names, output_path, title, ylabel)


def plot_channel_feature_heatmap(matrix, channel_labels, feature_names, output_path, title):
    fig_width = max(12.0, 0.18 * len(feature_names))
    fig_height = max(6.0, 0.45 * len(channel_labels))
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    im = ax.imshow(matrix, aspect="auto", origin="lower", cmap="magma")
    ax.set_xticks(np.arange(len(feature_names)))
    ax.set_xticklabels(feature_names, rotation=90, fontsize=7)
    ax.set_yticks(np.arange(len(channel_labels)))
    ax.set_yticklabels(channel_labels, fontsize=9)
    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax, shrink=0.85)
    cbar.ax.tick_params(labelsize=8)
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def plot_matrix(matrix, feature_names, output_path, title, vmin=None, vmax=None, cmap="coolwarm"):
    fig_size = max(10, int(len(feature_names) * 0.23))
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))
    im = ax.imshow(matrix, aspect="auto", origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks(np.arange(len(feature_names)))
    ax.set_xticklabels(feature_names, rotation=90, fontsize=7)
    ax.set_yticks(np.arange(len(feature_names)))
    ax.set_yticklabels(feature_names, fontsize=7)
    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax, shrink=0.82)
    cbar.ax.tick_params(labelsize=8)
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def active_feature_metadata(active_feature_groups):
    active_group_set = set(active_feature_groups)
    active_indices = []
    for group_name, feature_indices in build_feature_groups():
        if group_name in active_group_set:
            active_indices.extend(feature_indices)
    active_indices = sorted(active_indices)
    active_names = [FEATURE_FIELD_NAMES[index] for index in active_indices]
    return active_indices, active_names


def select_feature_values(values, active_feature_indices):
    return np.asarray(values, dtype=np.float64)[active_feature_indices]


def save_gradient_analysis(
    model,
    features,
    labels,
    weights,
    batch_size,
    device,
    output_dir,
    class_names,
    gradient_views,
    active_feature_indices,
    active_feature_names,
):
    gradient_dir = output_dir / "gradients"
    gradient_dir.mkdir(parents=True, exist_ok=True)

    by_target = []
    for class_index, class_name in enumerate(class_names):
        target_imp = compute_grad_importance(
            model,
            features,
            weights,
            batch_size,
            device,
            target_index=class_index,
        )
        by_target.append((class_name, select_feature_values(target_imp, active_feature_indices)))
    plot_feature_bar(
        by_target,
        active_feature_names,
        gradient_dir / "gradient_abs_by_target.png",
        "Absolute Gradient Importance by Target Score (d target score / d input)",
        "weighted mean |d score / d x|",
    )

    overall = None
    binary_margin = None
    if gradient_views == "full":
        overall = compute_grad_importance(model, features, weights, batch_size, device, target_index=None)
        plot_feature_bar(
            [("overall", select_feature_values(overall, active_feature_indices))],
            active_feature_names,
            gradient_dir / "gradient_abs_overall.png",
            "Absolute Gradient Importance on Predicted Targets",
            "weighted mean |d score / d x|",
        )

        by_truth = []
        for class_index, class_name in enumerate(class_names):
            mask = labels == class_index
            if not np.any(mask):
                continue
            truth_imp = compute_grad_importance(
                model,
                features[mask],
                weights[mask],
                batch_size,
                device,
                target_index=class_index,
            )
            by_truth.append((class_name, select_feature_values(truth_imp, active_feature_indices)))
        plot_feature_bar(
            by_truth,
            active_feature_names,
            gradient_dir / "gradient_abs_true_class.png",
            "Absolute Gradient Importance on True-Class Subsets",
            "weighted mean |d score / d x|",
        )

    if len(class_names) == 2:
        binary_margin = compute_binary_margin_importance(model, features, weights, batch_size, device)
        plot_feature_bar(
            [("binary_margin", select_feature_values(binary_margin, active_feature_indices))],
            active_feature_names,
            gradient_dir / "gradient_abs_binary_margin.png",
            f"Absolute Gradient Importance on logit({class_names[1]}) - logit({class_names[0]})",
            "weighted mean |d margin / d x|",
        )

    csv_path = gradient_dir / "gradient_importance_summary.csv"
    header = ["feature"]
    if overall is not None:
        header.append("overall")
    if binary_margin is not None:
        header.append("binary_margin")
    header.extend(f"target_{name}" for name in class_names)
    target_map = {label: values for label, values in by_target}
    with open(csv_path, "w", encoding="utf-8") as handle:
        handle.write(",".join(header) + "\n")
        for position, (active_index, feature_name) in enumerate(zip(active_feature_indices, active_feature_names)):
            row = [feature_name]
            if overall is not None:
                row.append(f"{overall[active_index]:.8g}")
            if binary_margin is not None:
                row.append(f"{binary_margin[active_index]:.8g}")
            row.extend(f"{target_map[name][position]:.8g}" for name in class_names)
            handle.write(",".join(row) + "\n")


def save_permutation_analysis(
    model,
    features,
    labels,
    weights,
    batch_size,
    output_dir,
    n_repeats,
    max_events,
    seed,
    active_feature_indices,
    active_feature_names,
):
    permutation_dir = output_dir / "permutation"
    permutation_dir.mkdir(parents=True, exist_ok=True)

    if features.shape[0] > max_events:
        rng = np.random.default_rng(seed)
        selected = rng.choice(features.shape[0], size=max_events, replace=False)
        subset_features = features[selected]
        subset_labels = labels[selected]
        subset_weights = weights[selected]
    else:
        subset_features = features
        subset_labels = labels
        subset_weights = weights

    device = next(model.parameters()).device
    importance, base_score = compute_permutation_importance(
        model,
        subset_features,
        subset_labels,
        subset_weights,
        batch_size,
        device,
        n_repeats=n_repeats,
        seed=seed,
        feature_indices=active_feature_indices,
    )
    plot_feature_bar(
        [("perm_drop", importance)],
        active_feature_names,
        permutation_dir / "permutation_importance.png",
        f"Permutation Importance (base weighted acc = {base_score:.4f})",
        "weighted accuracy drop",
    )
    with open(permutation_dir / "permutation_importance_summary.csv", "w", encoding="utf-8") as handle:
        handle.write("feature,accuracy_drop\n")
        for feature_name, value in zip(active_feature_names, importance):
            handle.write(f"{feature_name},{value:.8g}\n")

    if len(np.unique(subset_labels)) == 2:
        by_class = []
        for class_index in sorted(np.unique(subset_labels)):
            class_mask = subset_labels == int(class_index)
            if not np.any(class_mask):
                continue
            class_importance, class_base_score = compute_permutation_importance(
                model,
                subset_features[class_mask],
                subset_labels[class_mask],
                subset_weights[class_mask],
                batch_size,
                device,
                n_repeats=n_repeats,
                seed=seed + int(class_index) + 1000,
                feature_indices=active_feature_indices,
            )
            by_class.append((int(class_index), class_importance, class_base_score))

        if by_class:
            plot_feature_bar(
                [(f"class_{class_index}", values) for class_index, values, _ in by_class],
                active_feature_names,
                permutation_dir / "permutation_importance_by_true_class.png",
                "Permutation Importance by True Class",
                "weighted accuracy drop",
            )
            with open(permutation_dir / "permutation_importance_by_true_class.csv", "w", encoding="utf-8") as handle:
                header = ["feature"] + [f"class_{class_index}" for class_index, _, _ in by_class]
                handle.write(",".join(header) + "\n")
                for feature_index, feature_name in enumerate(active_feature_names):
                    row = [feature_name]
                    for _, values, _ in by_class:
                        row.append(f"{values[feature_index]:.8g}")
                    handle.write(",".join(row) + "\n")


def save_group_permutation_analysis(
    model,
    features,
    labels,
    weights,
    batch_size,
    output_dir,
    n_repeats,
    max_events,
    seed,
    active_feature_groups,
    active_feature_indices,
    active_feature_names,
    feature_layout,
):
    permutation_dir = output_dir / "group_permutation"
    permutation_dir.mkdir(parents=True, exist_ok=True)

    if features.shape[0] > max_events:
        rng = np.random.default_rng(seed)
        selected = rng.choice(features.shape[0], size=max_events, replace=False)
        subset_features = features[selected]
        subset_labels = labels[selected]
        subset_weights = weights[selected]
    else:
        subset_features = features
        subset_labels = labels
        subset_weights = weights

    active_group_set = set(active_feature_groups)
    groups = []
    compact_index_map = {feature_index: position for position, feature_index in enumerate(active_feature_indices)}
    for group_name, feature_indices in build_feature_groups():
        if group_name not in active_group_set:
            continue
        if feature_layout == "masked":
            groups.append((group_name, list(feature_indices)))
        else:
            groups.append((group_name, [compact_index_map[index] for index in feature_indices]))
    group_names = [group_name for group_name, _ in groups]
    device = next(model.parameters()).device
    importance, base_score = compute_group_permutation_importance(
        model,
        subset_features,
        subset_labels,
        subset_weights,
        batch_size,
        device,
        groups,
        n_repeats=n_repeats,
        seed=seed,
    )

    plot_named_bar(
        [("group_perm_drop", importance)],
        group_names,
        permutation_dir / "group_permutation_importance.png",
        f"Grouped Permutation Importance (base weighted acc = {base_score:.4f})",
        "weighted accuracy drop",
    )

    with open(permutation_dir / "group_permutation_importance_summary.csv", "w", encoding="utf-8") as handle:
        handle.write("group_name,n_features,features,accuracy_drop\n")
        for (group_name, feature_indices), value in zip(groups, importance):
            if feature_layout == "masked":
                feature_list = "|".join(FEATURE_FIELD_NAMES[index] for index in feature_indices)
            else:
                feature_list = "|".join(active_feature_names[index] for index in feature_indices)
            handle.write(f"{group_name},{len(feature_indices)},{feature_list},{value:.8g}\n")

    if len(np.unique(subset_labels)) == 2:
        by_class = []
        for class_index in sorted(np.unique(subset_labels)):
            class_mask = subset_labels == int(class_index)
            if not np.any(class_mask):
                continue
            class_importance, class_base_score = compute_group_permutation_importance(
                model,
                subset_features[class_mask],
                subset_labels[class_mask],
                subset_weights[class_mask],
                batch_size,
                device,
                groups,
                n_repeats=n_repeats,
                seed=seed + int(class_index) + 2000,
            )
            by_class.append((int(class_index), class_importance, class_base_score))

        if by_class:
            plot_named_bar(
                [(f"class_{class_index}", values) for class_index, values, _ in by_class],
                group_names,
                permutation_dir / "group_permutation_importance_by_true_class.png",
                "Grouped Permutation Importance by True Class",
                "weighted accuracy drop",
            )
            with open(
                permutation_dir / "group_permutation_importance_by_true_class.csv",
                "w",
                encoding="utf-8",
            ) as handle:
                header = ["group_name"] + [f"class_{class_index}" for class_index, _, _ in by_class]
                handle.write(",".join(header) + "\n")
                for group_index, group_name in enumerate(group_names):
                    row = [group_name]
                    for _, values, _ in by_class:
                        row.append(f"{values[group_index]:.8g}")
                    handle.write(",".join(row) + "\n")


def save_channel_permutation_analysis(
    model,
    features,
    labels,
    channel_ids,
    weights,
    batch_size,
    output_dir,
    n_repeats,
    max_events,
    min_events,
    seed,
    channel_id_to_info,
    active_feature_indices,
    active_feature_names,
):
    if channel_ids.size == 0 or not channel_id_to_info:
        return

    permutation_dir = output_dir / "channel_permutation"
    permutation_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    channel_rows = []

    for channel_id in sorted(np.unique(channel_ids)):
        channel_info = channel_id_to_info.get(int(channel_id))
        if channel_info is None:
            continue
        channel_name = channel_info["channel_name"]
        category_name = channel_info["class_name"]
        mask = channel_ids == int(channel_id)
        n_events = int(np.sum(mask))
        if n_events < min_events:
            continue

        subset_features = features[mask]
        subset_labels = labels[mask]
        subset_weights = weights[mask]
        if subset_features.shape[0] > max_events:
            selected = rng.choice(subset_features.shape[0], size=max_events, replace=False)
            subset_features = subset_features[selected]
            subset_labels = subset_labels[selected]
            subset_weights = subset_weights[selected]

        importance, base_score = compute_permutation_importance(
            model,
            subset_features,
            subset_labels,
            subset_weights,
            batch_size,
            next(model.parameters()).device,
            n_repeats=n_repeats,
            seed=seed + int(channel_id),
            feature_indices=active_feature_indices,
        )
        plot_feature_bar(
            [(channel_name, importance)],
            active_feature_names,
            permutation_dir / f"{channel_name}_permutation_importance.png",
            f"{channel_name} permutation importance (base weighted acc = {base_score:.4f})",
            "weighted accuracy drop",
        )
        with open(permutation_dir / f"{channel_name}_permutation_importance.csv", "w", encoding="utf-8") as handle:
            handle.write("feature,accuracy_drop\n")
            for feature_name, value in zip(active_feature_names, importance):
                handle.write(f"{feature_name},{value:.8g}\n")

        total = np.sum(importance, dtype=np.float64)
        normalized = importance / total if total > 0 else importance
        channel_rows.append(
            {
                "channel_name": channel_name,
                "category_name": category_name,
                "n_events": subset_features.shape[0],
                "base_score": base_score,
                "importance": importance,
                "normalized": normalized,
            }
        )

    if not channel_rows:
        return

    ordered_rows = sorted(channel_rows, key=lambda row: (row["category_name"], row["channel_name"]))
    heatmap = np.stack([row["normalized"] for row in ordered_rows], axis=0)
    labels_for_plot = [
        f"{row['channel_name']} ({row['category_name']}, n={row['n_events']})"
        for row in ordered_rows
    ]
    plot_channel_feature_heatmap(
        heatmap,
        labels_for_plot,
        active_feature_names,
        permutation_dir / "channel_permutation_heatmap_normalized.png",
        "Channel-wise permutation importance (row-normalized)",
    )

    with open(permutation_dir / "channel_permutation_summary.csv", "w", encoding="utf-8") as handle:
        header = ["channel", "category", "n_events", "base_weighted_acc", *active_feature_names]
        handle.write(",".join(header) + "\n")
        for row in ordered_rows:
            values = ",".join(f"{value:.8g}" for value in row["importance"])
            handle.write(
                f"{row['channel_name']},{row['category_name']},{row['n_events']},{row['base_score']:.8g},{values}\n"
            )


def save_covariance_analysis(features, labels, weights, output_dir, class_names, active_feature_names):
    covariance_dir = output_dir / "covariance"
    covariance_dir.mkdir(parents=True, exist_ok=True)

    overall_cov = weighted_covariance(features, weights)
    overall_corr = covariance_to_correlation(overall_cov)
    plot_matrix(
        overall_cov,
        active_feature_names,
        covariance_dir / "covariance_overall.png",
        "Weighted Feature Covariance",
        cmap="viridis",
    )
    plot_matrix(
        overall_corr,
        active_feature_names,
        covariance_dir / "correlation_overall.png",
        "Weighted Feature Correlation",
        vmin=-1.0,
        vmax=1.0,
        cmap="coolwarm",
    )

    for class_index, class_name in enumerate(class_names):
        mask = labels == class_index
        if not np.any(mask):
            continue
        class_cov = weighted_covariance(features[mask], weights[mask])
        class_corr = covariance_to_correlation(class_cov)
        plot_matrix(
            class_cov,
            active_feature_names,
            covariance_dir / f"covariance_{class_name}.png",
            f"Weighted Feature Covariance: {class_name}",
            cmap="viridis",
        )
        plot_matrix(
            class_corr,
            active_feature_names,
            covariance_dir / f"correlation_{class_name}.png",
            f"Weighted Feature Correlation: {class_name}",
            vmin=-1.0,
            vmax=1.0,
            cmap="coolwarm",
        )
        plot_matrix(
            class_corr - overall_corr,
            active_feature_names,
            covariance_dir / f"correlation_diff_{class_name}_minus_overall.png",
            f"Correlation Difference: {class_name} - overall",
            vmin=-1.0,
            vmax=1.0,
            cmap="coolwarm",
        )


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    report_path, dataset_report = load_dataset_contract(args.dataset, args.dataset_report)
    validate_mixed_dataset_file(args.dataset, dataset_report)
    source_label_names = list(dataset_report["label_names"])
    requested_subset = resolve_label_subset(args.label_subset, source_label_names)
    channel_id_to_info = load_channel_id_to_info(dataset_report)
    device = resolve_device(args.device)
    base_contract = build_base_data_contract(dataset_report)
    model_state_dict, saved_contract = load_model_checkpoint(args.model, device)
    effective_subset = saved_contract["task_label_indices"] if requested_subset is None else requested_subset
    active_feature_groups = list(saved_contract["active_feature_groups"])
    task_contract = build_task_data_contract(
        dataset_report,
        effective_subset,
        feature_groups=active_feature_groups,
    )
    task_contract["feature_layout"] = saved_contract["feature_layout"]
    validate_checkpoint_contract(saved_contract, task_contract, args.model)
    class_names = list(task_contract["task_label_names"])
    active_feature_indices, active_feature_names = active_feature_metadata(active_feature_groups)
    feature_layout = saved_contract["feature_layout"]
    features, raw_weights, labels, is_same_sign, channel_ids = loadDS(args.dataset)
    features, raw_weights, labels, is_same_sign, channel_ids = filter_dataset_for_task(
        features, raw_weights, labels, is_same_sign, channel_ids, effective_subset
    )
    features = apply_feature_mask(features, active_feature_groups)
    active_features = select_active_features(features, active_feature_groups)
    weights = resolve_analysis_weights(raw_weights, labels, args.weight_mode)
    model_inputs = features if feature_layout == "masked" else active_features
    model = SimpleMLP(input_dim=model_inputs.shape[1], n_classes=len(class_names)).to(device)
    model.load_state_dict(model_state_dict)
    model.eval()

    save_gradient_analysis(
        model=model,
        features=model_inputs,
        labels=labels,
        weights=weights,
        batch_size=args.batch_size,
        device=device,
        output_dir=output_dir,
        class_names=class_names,
        gradient_views=args.gradient_views,
        active_feature_indices=active_feature_indices if feature_layout == "masked" else list(range(len(active_feature_names))),
        active_feature_names=active_feature_names,
    )
    save_permutation_analysis(
        model=model,
        features=model_inputs,
        labels=labels,
        weights=weights,
        batch_size=args.batch_size,
        output_dir=output_dir,
        n_repeats=args.perm_repeats,
        max_events=args.perm_max_events,
        seed=args.seed,
        active_feature_indices=active_feature_indices if feature_layout == "masked" else list(range(len(active_feature_names))),
        active_feature_names=active_feature_names,
    )
    save_group_permutation_analysis(
        model=model,
        features=model_inputs,
        labels=labels,
        weights=weights,
        batch_size=args.batch_size,
        output_dir=output_dir,
        n_repeats=args.perm_repeats,
        max_events=args.perm_max_events,
        seed=args.seed,
        active_feature_groups=active_feature_groups,
        active_feature_indices=active_feature_indices,
        active_feature_names=active_feature_names,
        feature_layout=feature_layout,
    )
    save_channel_permutation_analysis(
        model=model,
        features=model_inputs,
        labels=labels,
        channel_ids=channel_ids,
        weights=weights,
        batch_size=args.batch_size,
        output_dir=output_dir,
        n_repeats=args.perm_repeats,
        max_events=args.perm_max_events,
        min_events=args.perm_min_events,
        seed=args.seed,
        channel_id_to_info=channel_id_to_info,
        active_feature_indices=active_feature_indices if feature_layout == "masked" else list(range(len(active_feature_names))),
        active_feature_names=active_feature_names,
    )
    save_covariance_analysis(
        features=active_features.astype(np.float64),
        labels=labels,
        weights=weights,
        output_dir=output_dir,
        class_names=class_names,
        active_feature_names=active_feature_names,
    )

    print(f"validated dataset report: {report_path}")
    print("analysis task labels:", ", ".join(class_names))
    print("active feature groups:", ", ".join(active_feature_groups))
    print("analysis outputs: gradients, permutation, group permutation, channel permutation, covariance")
    print(f"saved output dir: {output_dir}")


if __name__ == "__main__":
    main()
