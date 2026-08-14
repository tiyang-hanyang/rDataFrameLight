import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


KNOWN_PREDICTION_FIELDS = {"rowIndex", "truth", "prediction", "weight", "isSameSign"}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare multiple simpleSignalModel run outputs: scores, gradients, and permutation studies."
    )
    parser.add_argument("--runs", nargs="+", required=True, help="Run directories containing model outputs.")
    parser.add_argument(
        "--labels",
        nargs="+",
        default=None,
        help="Optional display labels for --runs. Defaults to directory names.",
    )
    parser.add_argument("--output-dir", required=True, help="Output directory for comparison products.")
    parser.add_argument(
        "--top-n",
        type=int,
        default=15,
        help="Number of top features to show in gradient/permutation comparison plots.",
    )
    return parser.parse_args()


def resolve_run_labels(run_paths, labels):
    if labels is None:
        return [path.name for path in run_paths]
    if len(labels) != len(run_paths):
        raise RuntimeError("--labels must have the same length as --runs.")
    return list(labels)


def load_predictions(run_dir):
    prediction_path = run_dir / "predictions.npy"
    if not prediction_path.is_file():
        raise RuntimeError(f"Missing prediction file: {prediction_path}")
    array = np.load(prediction_path, allow_pickle=False)
    score_fields = [field for field in (array.dtype.names or ()) if field not in KNOWN_PREDICTION_FIELDS]
    if len(score_fields) < 2:
        raise RuntimeError(f"{prediction_path} does not contain at least two score fields.")
    if "weight" not in (array.dtype.names or ()):
        raise RuntimeError(f"{prediction_path} is missing the weight field.")
    return array, score_fields


def load_csv_table(csv_path, key_field):
    if not csv_path.is_file():
        raise RuntimeError(f"Missing CSV file: {csv_path}")
    with open(csv_path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    if not rows:
        return {}, []
    columns = [name for name in reader.fieldnames if name != key_field]
    table = {row[key_field]: row for row in rows}
    return table, columns


def load_optional_csv_table(csv_path, key_field):
    if not csv_path.is_file():
        return None, []
    return load_csv_table(csv_path, key_field)


def weighted_accuracy(predictions):
    truth = predictions["truth"].astype(np.int64)
    pred = predictions["prediction"].astype(np.int64)
    weights = np.abs(predictions["weight"].astype(np.float64))
    weight_sum = np.sum(weights, dtype=np.float64)
    if weight_sum <= 0:
        return float(np.mean(pred == truth)) if truth.size > 0 else 0.0
    return float(np.sum(weights * (pred == truth), dtype=np.float64) / weight_sum)


def unweighted_accuracy(predictions):
    truth = predictions["truth"].astype(np.int64)
    pred = predictions["prediction"].astype(np.int64)
    return float(np.mean(pred == truth)) if truth.size > 0 else 0.0


def normalize_signed_class_weights(weights):
    weight_sum = np.sum(weights, dtype=np.float64)
    if np.abs(weight_sum) > 0:
        return weights / weight_sum, weight_sum
    return weights, 0.0


def format_class_label(name):
    mapping = {
        "TTHH": r"$t\bar{t}HH$",
        "ttbar_ttbb": r"$t\bar{t}(b\bar{b})$",
        "tt_b": r"$t\bar{t}(b\bar{b})$",
        "ttX_like": r"$t\bar{t}+X$",
        "ttW_like": r"$t\bar{t}+X$",
        "TTTT": r"$t\bar{t}t\bar{t}$",
    }
    return mapping.get(name, name)


def weighted_hist_with_error(values, weights, bins):
    values = np.asarray(values, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    clipped = values.copy()
    low_edge = np.nextafter(bins[0], np.inf)
    high_edge = np.nextafter(bins[-1], -np.inf)
    clipped[clipped < bins[0]] = low_edge
    clipped[clipped >= bins[-1]] = high_edge
    hist = np.histogram(clipped, bins=bins, weights=weights)[0].astype(np.float64)
    err2 = np.histogram(clipped, bins=bins, weights=np.square(weights))[0].astype(np.float64)
    return hist, np.sqrt(np.maximum(err2, 0.0))


def plot_run_score_shapes(output_dir, run_payloads, score_fields):
    score_dir = output_dir / "scores"
    score_dir.mkdir(parents=True, exist_ok=True)
    bins = np.linspace(0.0, 1.0, 31)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])

    for score_field in score_fields:
        for truth_index, truth_name in enumerate(score_fields):
            fig, ax = plt.subplots(figsize=(8, 6))
            plotted = False
            for payload in run_payloads:
                predictions = payload["predictions"]
                mask = predictions["truth"].astype(np.int64) == truth_index
                if not np.any(mask):
                    continue
                weights, weight_sum = normalize_signed_class_weights(
                    predictions["weight"][mask].astype(np.float64)
                )
                if np.abs(weight_sum) > 0:
                    hist, err = weighted_hist_with_error(predictions[score_field][mask], weights, bins)
                else:
                    unit_weights = np.ones(np.count_nonzero(mask), dtype=np.float64)
                    hist, err = weighted_hist_with_error(predictions[score_field][mask], unit_weights, bins)
                stairs = ax.stairs(
                    hist,
                    bins,
                    linewidth=1.8,
                    label=payload["label"],
                )
                ax.errorbar(
                    bin_centers,
                    hist,
                    yerr=err,
                    fmt="o",
                    color=stairs.get_edgecolor(),
                    markerfacecolor="white",
                    markeredgecolor=stairs.get_edgecolor(),
                    markersize=2.8,
                    elinewidth=1.0,
                    capsize=2.2,
                    zorder=3,
                )
                plotted = True
            ax.set_xlabel(f"{format_class_label(score_field)} score")
            ax.set_ylabel("Normalized weighted entries")
            ax.text(
                0.03,
                0.97,
                f"{format_class_label(truth_name)} test dataset",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=13,
                fontweight="bold",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 3.0},
            )
            ax.set_xlim(float(bins[0]), float(bins[-1]))
            ax.grid(alpha=0.25)
            if plotted:
                ax.legend(fontsize=9)
            fig.tight_layout()
            fig.savefig(score_dir / f"{truth_name}_truth__{score_field}_score_runs.png")
            plt.close(fig)


def save_score_summary(output_path, run_payloads, score_fields):
    with open(output_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["run", "truth_class", "score_field", "weighted_mean", "unweighted_mean"],
        )
        writer.writeheader()
        for payload in run_payloads:
            predictions = payload["predictions"]
            truth = predictions["truth"].astype(np.int64)
            for truth_index, truth_name in enumerate(score_fields):
                mask = truth == truth_index
                if not np.any(mask):
                    continue
                weights = predictions["weight"][mask].astype(np.float64)
                weight_sum = np.sum(weights, dtype=np.float64)
                for score_field in score_fields:
                    score_values = predictions[score_field][mask].astype(np.float64)
                    weighted_mean = (
                        float(np.sum(weights * score_values, dtype=np.float64) / weight_sum)
                        if np.abs(weight_sum) > 0
                        else float(np.mean(score_values))
                    )
                    writer.writerow(
                        {
                            "run": payload["label"],
                            "truth_class": truth_name,
                            "score_field": score_field,
                            "weighted_mean": f"{weighted_mean:.8g}",
                            "unweighted_mean": f"{float(np.mean(score_values)):.8g}",
                        }
                    )


def plot_metric_comparison(output_path, metric_label, run_payloads, table_key, value_key, top_n):
    value_map = {}
    for payload in run_payloads:
        table = payload[table_key]
        for item_name, row in table.items():
            if value_key not in row or row[value_key] == "":
                continue
            value_map.setdefault(item_name, {})[payload["label"]] = float(row[value_key])

    if not value_map:
        return

    ranking = sorted(
        value_map.items(),
        key=lambda item: np.mean(list(item[1].values()), dtype=np.float64),
        reverse=True,
    )[:top_n]
    selected_names = [name for name, _ in ranking]

    fig_height = max(6.0, 0.35 * len(selected_names))
    fig, ax = plt.subplots(figsize=(11, fig_height))
    y_pos = np.arange(len(selected_names))
    n_runs = len(run_payloads)
    bar_height = 0.8 / max(n_runs, 1)

    for run_index, payload in enumerate(run_payloads):
        values = [value_map.get(name, {}).get(payload["label"], 0.0) for name in selected_names]
        offset = (run_index - (n_runs - 1) / 2.0) * bar_height
        ax.barh(y_pos + offset, values, height=bar_height, label=payload["label"], alpha=0.85)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(selected_names, fontsize=8)
    ax.set_xlabel(metric_label)
    ax.grid(alpha=0.25, axis="x")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def save_wide_summary(output_path, row_key_name, run_payloads, table_key, value_key):
    rows = {}
    for payload in run_payloads:
        table = payload[table_key]
        for item_name, row in table.items():
            if value_key not in row or row[value_key] == "":
                continue
            rows.setdefault(item_name, {})[payload["label"]] = float(row[value_key])

    if not rows:
        return

    with open(output_path, "w", encoding="utf-8", newline="") as handle:
        fieldnames = [row_key_name] + [payload["label"] for payload in run_payloads]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for item_name in sorted(rows):
            record = {row_key_name: item_name}
            record.update(rows[item_name])
            writer.writerow(record)


def auto_gradient_column(columns):
    for preferred in ("binary_margin", "overall"):
        if preferred in columns:
            return preferred
    target_columns = [column for column in columns if column.startswith("target_")]
    if target_columns:
        return target_columns[-1]
    if columns:
        return columns[0]
    raise RuntimeError("Gradient summary does not contain any metric columns.")


def main():
    args = parse_args()
    run_dirs = [Path(run).resolve() for run in args.runs]
    run_labels = resolve_run_labels(run_dirs, args.labels)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    run_payloads = []
    reference_score_fields = None
    reference_gradient_column = None
    notices = []

    for run_dir, run_label in zip(run_dirs, run_labels):
        predictions, score_fields = load_predictions(run_dir)
        if reference_score_fields is None:
            reference_score_fields = score_fields
        elif score_fields != reference_score_fields:
            raise RuntimeError(
                f"{run_dir}: score fields {score_fields} do not match reference {reference_score_fields}."
            )

        gradient_table, gradient_columns = load_optional_csv_table(
            run_dir / "analysis" / "gradients" / "gradient_importance_summary.csv",
            "feature",
        )
        gradient_column = None
        if gradient_table is not None:
            gradient_column = auto_gradient_column(gradient_columns)
            if reference_gradient_column is None:
                reference_gradient_column = gradient_column
            elif gradient_column != reference_gradient_column:
                raise RuntimeError(
                    f"{run_dir}: selected gradient column '{gradient_column}' does not match reference "
                    f"'{reference_gradient_column}'."
                )
        else:
            notices.append(f"[{run_label}] no gradient summary; score comparison only for this run")
            gradient_table = {}

        permutation_table, _ = load_optional_csv_table(
            run_dir / "analysis" / "permutation" / "permutation_importance_summary.csv",
            "feature",
        )
        if permutation_table is None:
            notices.append(f"[{run_label}] no permutation summary; skipping permutation comparison for this run")
            permutation_table = {}

        group_permutation_table, _ = load_optional_csv_table(
            run_dir / "analysis" / "group_permutation" / "group_permutation_importance_summary.csv",
            "group_name",
        )
        if group_permutation_table is None:
            notices.append(
                f"[{run_label}] no group permutation summary; skipping grouped permutation comparison for this run"
            )
            group_permutation_table = {}

        run_payloads.append(
            {
                "label": run_label,
                "run_dir": run_dir,
                "predictions": predictions,
                "gradient_table": gradient_table,
                "permutation_table": permutation_table,
                "group_permutation_table": group_permutation_table,
            }
        )

    summary_csv = output_dir / "run_summary.csv"
    with open(summary_csv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["run", "unweighted_accuracy", "weighted_accuracy_abs"],
        )
        writer.writeheader()
        for payload in run_payloads:
            writer.writerow(
                {
                    "run": payload["label"],
                    "unweighted_accuracy": f"{unweighted_accuracy(payload['predictions']):.8g}",
                    "weighted_accuracy_abs": f"{weighted_accuracy(payload['predictions']):.8g}",
                }
            )

    plot_run_score_shapes(output_dir, run_payloads, reference_score_fields)
    save_score_summary(output_dir / "score_summary.csv", run_payloads, reference_score_fields)
    if reference_gradient_column is not None:
        save_wide_summary(
            output_dir / "gradient_comparison.csv",
            "feature",
            run_payloads,
            "gradient_table",
            reference_gradient_column,
        )
        plot_metric_comparison(
            output_dir / "gradient_top_features.png",
            f"gradient ({reference_gradient_column})",
            run_payloads,
            "gradient_table",
            reference_gradient_column,
            args.top_n,
        )

    save_wide_summary(
        output_dir / "permutation_comparison.csv",
        "feature",
        run_payloads,
        "permutation_table",
        "accuracy_drop",
    )
    save_wide_summary(
        output_dir / "group_permutation_comparison.csv",
        "group_name",
        run_payloads,
        "group_permutation_table",
        "accuracy_drop",
    )

    plot_metric_comparison(
        output_dir / "permutation_top_features.png",
        "permutation accuracy drop",
        run_payloads,
        "permutation_table",
        "accuracy_drop",
        args.top_n,
    )
    plot_metric_comparison(
        output_dir / "group_permutation_all_groups.png",
        "group permutation accuracy drop",
        run_payloads,
        "group_permutation_table",
        "accuracy_drop",
        20,
    )

    print("compared runs:", ", ".join(payload["label"] for payload in run_payloads))
    print("score fields:", ", ".join(reference_score_fields))
    print("gradient column:", reference_gradient_column if reference_gradient_column is not None else "N/A")
    for notice in notices:
        print("note:", notice)
    print(f"saved output dir: {output_dir}")


if __name__ == "__main__":
    main()
