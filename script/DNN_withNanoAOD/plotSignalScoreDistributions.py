import argparse
import json
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt

from prepareSignalDatasets import (
    CHANNEL_INDEX,
    DEFAULT_MIX_REPORT_NAME,
    load_json,
    resolve_dataset_report_path,
    structured_to_matrix,
    validate_dataset_report,
    validate_mixed_dataset_file,
)
from simpleSignalModel import build_task_data_contract, resolve_label_subset
from signal_class_config import (
    DEFAULT_SCHEMA_NAME,
    get_category_color,
    get_channel_color,
    get_signal_schema,
    schema_help_text,
)
BIN_EDGES = np.asarray([0.0, 0.25, 0.5, 0.75, 1.0], dtype=np.float32)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot score distributions by category and channel."
    )
    parser.add_argument(
        "--inputs",
        nargs="+",
        default=None,
        help="Per-channel structured prediction .npy files.",
    )
    parser.add_argument(
        "--dataset",
        default=None,
        help="Dataset .npy used for testing, required together with --predictions for mixed samples.",
    )
    parser.add_argument(
        "--predictions",
        default=None,
        help="Structured prediction .npy saved by simpleSignalModel.py test --save-predictions.",
    )
    parser.add_argument(
        "--channel-map",
        default=None,
        help="Optional channel_map.json for decoding channelId in --dataset mode.",
    )
    parser.add_argument(
        "--dataset-report",
        default=None,
        help=(
            f"Optional {DEFAULT_MIX_REPORT_NAME} from prepareSignalDatasets.py. "
            "If omitted, it is resolved next to the dataset file."
        ),
    )
    parser.add_argument(
        "--label-subset",
        nargs="+",
        default=None,
        help="Optional sub-task labels, e.g. '--label-subset 0 1' or '--label-subset TTHH tt_b'.",
    )
    parser.add_argument(
        "--output-dir",
        default="signal_score_plots",
        help="Directory to store output figures.",
    )
    parser.add_argument(
        "--class-schema",
        default=DEFAULT_SCHEMA_NAME,
        help=schema_help_text(),
    )
    return parser.parse_args()


def load_prediction(prediction_path, score_fields):
    array = np.load(prediction_path, allow_pickle=False)
    missing = [field for field in score_fields if field not in (array.dtype.names or ())]
    if missing:
        raise RuntimeError(
            f"{prediction_path} is missing score fields: " + ", ".join(missing)
        )
    return array


def load_dataset(dataset_path):
    matrix = structured_to_matrix(np.load(dataset_path, allow_pickle=False))
    if matrix.shape[1] <= CHANNEL_INDEX:
        raise RuntimeError(f"{dataset_path} is missing channelId information.")
    return matrix


def load_dataset_contract(dataset_path, dataset_report_path=None):
    report_path = resolve_dataset_report_path(dataset_path, explicit_report_path=dataset_report_path)
    report = load_json(report_path)
    validate_dataset_report(report, str(report_path))
    return report_path, report


def infer_task_names_from_prediction_fields(prediction_array, source_label_names):
    score_fields = [field for field in (prediction_array.dtype.names or ()) if field in source_label_names]
    if len(score_fields) < 2:
        raise RuntimeError("Unable to infer plotting task labels from prediction score fields.")
    return score_fields


def default_channel_map(schema):
    return {channel.channel_name: index for index, channel in enumerate(schema.channels)}


def invert_channel_map(channel_map):
    return {int(channel_index): channel_name for channel_name, channel_index in channel_map.items()}


def load_channel_id_to_name(channel_map_path, schema):
    if channel_map_path is None:
        return invert_channel_map(default_channel_map(schema))

    with open(channel_map_path, "r", encoding="utf-8") as handle:
        raw_map = json.load(handle)
    normalized_map = {str(channel_name): int(channel_index) for channel_name, channel_index in raw_map.items()}
    return invert_channel_map(normalized_map)


def load_channel_id_to_name_from_report(channel_map_path, schema, dataset_report):
    if channel_map_path is not None:
        return load_channel_id_to_name(channel_map_path, schema)
    raw_info = dataset_report["labeling"]["channel_info"]
    return {
        int(channel_id): info["channel_name"]
        for channel_id, info in raw_info.items()
    }


def get_weights(array):
    if array.dtype.names is not None and "weight" in array.dtype.names:
        return array["weight"].astype(np.float64)
    return np.ones(array.shape[0], dtype=np.float64)


def clip_to_bins(values, bins):
    clipped = np.asarray(values, dtype=np.float64).copy()
    low_edge = np.nextafter(bins[0], np.inf)
    high_edge = np.nextafter(bins[-1], -np.inf)
    clipped[clipped < bins[0]] = low_edge
    clipped[clipped >= bins[-1]] = high_edge
    return clipped


def weighted_hist(scores, weights):
    clipped_scores = clip_to_bins(scores, BIN_EDGES)
    hist = np.histogram(clipped_scores, bins=BIN_EDGES, weights=weights)[0].astype(np.float64)
    err2 = np.histogram(clipped_scores, bins=BIN_EDGES, weights=np.square(weights))[0].astype(np.float64)
    return hist, np.sqrt(np.maximum(err2, 0.0))


def normalize_hist_and_error(hist, err, mode):
    normalized = hist.astype(np.float64).copy()
    normalized_err = err.astype(np.float64).copy()
    if mode == "none":
        return normalized, normalized_err
    if mode == "per_hist":
        scale = np.sum(normalized, dtype=np.float64)
        if np.abs(scale) > 0:
            normalized /= scale
            normalized_err /= scale
        return normalized, normalized_err
    raise RuntimeError(f"Unknown normalization mode: {mode}")


def make_series(label, scores, weights, color):
    return {
        "label": label,
        "scores": np.asarray(scores, dtype=np.float64),
        "weights": np.asarray(weights, dtype=np.float64),
        "color": color,
    }


def plot_series(output_path, title, xlabel, ylabel, series_list, normalization):
    fig, ax = plt.subplots(figsize=(8, 6))
    bin_centers = 0.5 * (BIN_EDGES[:-1] + BIN_EDGES[1:])
    histograms = []
    for series in series_list:
        hist, err = weighted_hist(series["scores"], series["weights"])
        histograms.append((series["label"], hist, err, series["color"]))

    if normalization == "global":
        total = np.sum([np.sum(hist, dtype=np.float64) for _, hist, _, _ in histograms], dtype=np.float64)
        global_scale = total if np.abs(total) > 0 else 1.0
    else:
        global_scale = None

    plotted = False
    for label, hist, err, color in histograms:
        if normalization == "per_hist":
            plotted_hist, plotted_err = normalize_hist_and_error(hist, err, "per_hist")
        elif normalization == "global":
            plotted_hist = hist / global_scale
            plotted_err = err / global_scale
        elif normalization == "none":
            plotted_hist = hist
            plotted_err = err
        else:
            raise RuntimeError(f"Unknown normalization mode: {normalization}")

        stairs = ax.stairs(plotted_hist, BIN_EDGES, linewidth=1.8, label=label, color=color)
        ax.errorbar(
            bin_centers,
            plotted_hist,
            yerr=plotted_err,
            fmt="none",
            color=stairs.get_edgecolor(),
            elinewidth=1.0,
            capsize=2,
        )
        plotted = True

    ax.set_xlim(float(BIN_EDGES[0]), float(BIN_EDGES[-1]))
    ax.set_xticks(BIN_EDGES)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(alpha=0.25)
    if plotted:
        ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def build_channel_payloads_from_dataset(dataset_path, prediction_path, channel_map_path, schema, class_names, dataset_report):
    dataset = load_dataset(dataset_path)
    predictions = load_prediction(prediction_path, class_names)
    if "rowIndex" not in (predictions.dtype.names or ()):
        raise RuntimeError(
            f"{prediction_path} is missing rowIndex. Re-run simpleSignalModel.py test --save-predictions with the updated script."
        )

    row_index = predictions["rowIndex"].astype(np.int64)
    if np.any(row_index < 0) or np.any(row_index >= dataset.shape[0]):
        raise RuntimeError("Prediction rowIndex is out of range for the provided dataset.")

    channel_ids = dataset[row_index, CHANNEL_INDEX].astype(np.int64)
    channel_id_to_name = load_channel_id_to_name_from_report(channel_map_path, schema, dataset_report)
    canonical_channel_id_to_name = {}
    for channel_id, mapped_name in channel_id_to_name.items():
        try:
            canonical_channel_id_to_name[channel_id] = schema.canonicalize_channel_name(mapped_name)
        except RuntimeError:
            continue

    channel_payloads = {}
    report_channel_info = {
        int(channel_id): info for channel_id, info in dataset_report["labeling"]["channel_info"].items()
    }
    seen_channels = []
    for channel_id in sorted(np.unique(channel_ids)):
        channel_info = report_channel_info.get(int(channel_id))
        if channel_info is None:
            continue
        channel_name = channel_info["channel_name"]
        category_name = channel_info["class_name"]
        if channel_name in seen_channels:
            continue
        seen_channels.append(channel_name)
        matching_ids = [
            channel_id
            for channel_id, canonical_name in canonical_channel_id_to_name.items()
            if canonical_name == channel_name
        ]
        if not matching_ids:
            continue
        mask = np.isin(channel_ids, matching_ids)
        if not np.any(mask):
            continue
        channel_payloads[channel_name] = {
            "category": category_name,
            "array": predictions[mask],
        }
    return channel_payloads


def plot_category_overall(output_dir, category_payloads, class_names):
    overall_dir = output_dir / "categories_overall"
    overall_dir.mkdir(parents=True, exist_ok=True)

    for score_field in class_names:
        series_list = []
        for category_name in class_names:
            arrays = [array for _, array in category_payloads[category_name]]
            if not arrays:
                continue
            merged = arrays[0] if len(arrays) == 1 else np.concatenate(arrays)
            series_list.append(
                make_series(
                    category_name,
                    merged[score_field],
                    get_weights(merged),
                    get_category_color(category_name),
                )
            )
        plot_series(
            overall_dir / f"{score_field}_categories.png",
            title=f"Category comparison in {score_field}",
            xlabel=f"{score_field} score",
            ylabel="Normalized weighted entries",
            series_list=series_list,
            normalization="per_hist",
        )


def plot_channel_score_shapes(output_dir, channel_payloads, class_names):
    channels_dir = output_dir / "channels_scores"
    channels_dir.mkdir(parents=True, exist_ok=True)

    for channel_name, payload in channel_payloads.items():
        array = payload["array"]
        weights = get_weights(array)
        series_list = [
            make_series(score_field, array[score_field], weights, get_category_color(score_field))
            for score_field in class_names
        ]
        plot_series(
            channels_dir / f"{channel_name}_scores.png",
            title=f"{channel_name} score shapes",
            xlabel="score",
            ylabel="Normalized weighted entries",
            series_list=series_list,
            normalization="per_hist",
        )


def plot_category_channel_breakdown(output_dir, category_payloads, class_names):
    breakdown_dir = output_dir / "category_channel_breakdown"
    global_dir = breakdown_dir / "global_normalized"
    per_hist_dir = breakdown_dir / "per_channel_normalized"
    global_dir.mkdir(parents=True, exist_ok=True)
    per_hist_dir.mkdir(parents=True, exist_ok=True)

    for category_name in class_names:
        members = category_payloads[category_name]
        if not members:
            continue
        for score_field in class_names:
            series_list = [
                make_series(channel_name, array[score_field], get_weights(array), get_channel_color(channel_name))
                for channel_name, array in members
            ]
            plot_series(
                global_dir / f"{category_name}_{score_field}_global.png",
                title=f"{category_name} channels in {score_field}",
                xlabel=f"{score_field} score",
                ylabel="Weighted entries / total category weight",
                series_list=series_list,
                normalization="global",
            )
            plot_series(
                per_hist_dir / f"{category_name}_{score_field}_per_channel.png",
                title=f"{category_name} channels in {score_field}",
                xlabel=f"{score_field} score",
                ylabel="Normalized weighted entries",
                series_list=series_list,
                normalization="per_hist",
            )


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    using_join_mode = args.dataset is not None or args.predictions is not None
    if using_join_mode:
        if args.dataset is None or args.predictions is None:
            raise RuntimeError("--dataset and --predictions must be provided together.")
        report_path, dataset_report = load_dataset_contract(args.dataset, args.dataset_report)
        validate_mixed_dataset_file(args.dataset, dataset_report)
        source_label_names = list(dataset_report["label_names"])
        requested_subset = resolve_label_subset(args.label_subset, source_label_names)
        raw_prediction = np.load(args.predictions, allow_pickle=False)
        if requested_subset is None:
            inferred_task_names = infer_task_names_from_prediction_fields(raw_prediction, source_label_names)
            effective_subset = [source_label_names.index(name) for name in inferred_task_names]
        else:
            effective_subset = requested_subset
        task_contract = build_task_data_contract(dataset_report, effective_subset)
        class_names = list(task_contract["task_label_names"])
        schema = get_signal_schema(args.class_schema)
        channel_payloads = build_channel_payloads_from_dataset(
            args.dataset,
            args.predictions,
            args.channel_map,
            schema,
            class_names,
            dataset_report,
        )
    else:
        if not args.inputs:
            raise RuntimeError("Provide either --inputs or both --dataset and --predictions.")
        schema = get_signal_schema(args.class_schema)
        class_names = None
        grouped = {}
        for input_path in args.inputs:
            raw_array = np.load(input_path, allow_pickle=False)
            if class_names is None:
                class_names = infer_task_names_from_prediction_fields(raw_array, list(schema.class_names))
            array = load_prediction(input_path, class_names)
            channel_name, category_name = schema.identify_channel(Path(input_path).stem)
            grouped.setdefault(channel_name, {"category": category_name, "arrays": []})
            grouped[channel_name]["arrays"].append(array)
        channel_payloads = {}
        for channel in schema.channels:
            channel_name = channel.channel_name
            if channel_name not in grouped:
                continue
            payload = grouped[channel_name]
            merged = payload["arrays"][0] if len(payload["arrays"]) == 1 else np.concatenate(payload["arrays"])
            channel_payloads[channel_name] = {
                "category": payload["category"],
                "array": merged,
            }

    if not channel_payloads:
        raise RuntimeError("No channels were found to plot.")

    category_payloads = {category_name: [] for category_name in class_names}
    channel_order = {channel.channel_name: index for index, channel in enumerate(schema.channels)}
    for channel_name, payload in sorted(channel_payloads.items(), key=lambda item: channel_order.get(item[0], 999)):
        if payload["category"] not in category_payloads:
            continue
        category_payloads[payload["category"]].append((channel_name, payload["array"]))
    plot_category_overall(output_dir, category_payloads, class_names)
    plot_channel_score_shapes(output_dir, channel_payloads, class_names)
    plot_category_channel_breakdown(output_dir, category_payloads, class_names)

    total_plots = len(class_names) + len(channel_payloads) + 2 * len(class_names) * len(class_names)
    if using_join_mode:
        print(f"validated dataset report: {report_path}")
        print("plot task labels:", ", ".join(class_names))
    print(f"saved output dir: {output_dir}")
    print(f"saved plots: {total_plots}")


if __name__ == "__main__":
    main()
