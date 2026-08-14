import argparse
import json
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt

from signalExtractionInputFromNanoAOD import FEATURE_FIELD_NAMES
from prepareSignalDatasets import CHANNEL_INDEX, FEATURE_DIM, LABEL_INDEX, WEIGHT_INDEX, structured_to_matrix
from signal_class_config import (
    DEFAULT_SCHEMA_NAME,
    get_category_color,
    get_channel_color,
    get_signal_schema,
    schema_help_text,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot input feature distributions by category and subchannel."
    )
    parser.add_argument("--dataset", required=True, help="Input dataset .npy, e.g. OS_train.npy.")
    parser.add_argument("--channel-map", required=True, help="channel_map.json from prepareSignalDatasets.py.")
    parser.add_argument("--output-dir", default="input_feature_plots", help="Output directory.")
    parser.add_argument("--nbins", type=int, default=20, help="Default number of bins for continuous features.")
    parser.add_argument(
        "--predictions",
        default=None,
        help="Optional structured predictions .npy from simpleSignalModel.py test --save-predictions.",
    )
    parser.add_argument(
        "--focus-channel",
        default=None,
        help="If set together with --predictions, compare feature distributions inside this channel by predicted class.",
    )
    parser.add_argument(
        "--group-by",
        choices=["prediction", "truth"],
        default="prediction",
        help="Grouping source used together with --predictions and --focus-channel.",
    )
    parser.add_argument(
        "--class-schema",
        default=DEFAULT_SCHEMA_NAME,
        help=schema_help_text(),
    )
    return parser.parse_args()


def load_channel_id_to_name(channel_map_path, schema):
    with open(channel_map_path, "r", encoding="utf-8") as handle:
        raw_map = json.load(handle)

    channel_id_to_name = {}
    for mapped_name, channel_id in raw_map.items():
        try:
            canonical_name = schema.canonicalize_channel_name(str(mapped_name))
        except RuntimeError:
            continue
        channel_id_to_name[int(channel_id)] = canonical_name
    return channel_id_to_name


def load_dataset(dataset_path):
    matrix = structured_to_matrix(np.load(dataset_path, allow_pickle=False))
    if matrix.shape[1] <= CHANNEL_INDEX:
        raise RuntimeError(f"{dataset_path} is missing channelId.")
    return matrix


def load_prediction_table(prediction_path, class_names):
    array = np.load(prediction_path, allow_pickle=False)
    required = {"rowIndex", "prediction", "truth", *class_names}
    missing = required.difference(array.dtype.names or ())
    if missing:
        raise RuntimeError(f"{prediction_path} is missing fields: {', '.join(sorted(missing))}")
    return array


def build_feature_bins(feature_name, nbins):
    if "dr" in feature_name:
        return np.linspace(0.0, 5.0, nbins + 1, dtype=np.float64)
    if "_m_" in feature_name or feature_name == "dimuon_mass":
        return np.linspace(0.0, 400.0, nbins + 1, dtype=np.float64)
    if feature_name == "jet_cent" or "btag" in feature_name:
        return np.linspace(0.0, 1.0, nbins + 1, dtype=np.float64)
    if feature_name.startswith("sublead_mu_p"):
        if feature_name.endswith("_px") or feature_name.endswith("_py"):
            return np.linspace(-200.0, 200.0, nbins + 1, dtype=np.float64)
        if feature_name.endswith("_pz"):
            return np.linspace(-400.0, 400.0, nbins + 1, dtype=np.float64)
    if feature_name.startswith("lead_mu_p") or feature_name.startswith("jet"):
        if feature_name.endswith("_px") or feature_name.endswith("_py"):
            return np.linspace(-400.0, 400.0, nbins + 1, dtype=np.float64)
        if feature_name.endswith("_pz"):
            return np.linspace(-500.0, 500.0, nbins + 1, dtype=np.float64)
        if feature_name.endswith("_mass"):
            return np.linspace(0.0, 50.0, nbins + 1, dtype=np.float64)
    if feature_name == "sumjet_pt":
        return np.linspace(0.0, 1500.0, nbins + 1, dtype=np.float64)
    if feature_name == "sumbjet_pt":
        return np.linspace(0.0, 1000.0, nbins + 1, dtype=np.float64)
    if feature_name == "nGoodJet":
        return np.arange(3.5, 10.5 + 1.0, 1.0, dtype=np.float64)
    if feature_name == "nBJet":
        return np.arange(2.5, 5.5 + 1.0, 1.0, dtype=np.float64)
    return np.linspace(-400.0, 400.0, nbins + 1, dtype=np.float64)


def clip_to_bins(values, bins):
    clipped = np.asarray(values, dtype=np.float64).copy()
    low_edge = np.nextafter(bins[0], np.inf)
    high_edge = np.nextafter(bins[-1], -np.inf)
    clipped[clipped < bins[0]] = low_edge
    clipped[clipped >= bins[-1]] = high_edge
    return clipped


def weighted_hist(values, weights, bins):
    clipped_values = clip_to_bins(values, bins)
    hist = np.histogram(clipped_values, bins=bins, weights=weights)[0].astype(np.float64)
    err2 = np.histogram(clipped_values, bins=bins, weights=np.square(weights))[0].astype(np.float64)
    return hist, np.sqrt(np.maximum(err2, 0.0))


def normalize_hist_and_error(hist, err, mode, global_scale=None):
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
    if mode == "global":
        scale = global_scale if global_scale is not None else 0.0
        if np.abs(scale) > 0:
            normalized /= scale
            normalized_err /= scale
        return normalized, normalized_err
    raise RuntimeError(f"Unknown normalization mode: {mode}")


def plot_hist_series(output_path, title, xlabel, ylabel, bins, series_list, normalization):
    fig, ax = plt.subplots(figsize=(8, 6))

    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    histograms = []
    for label, values, weights, color in series_list:
        hist, err = weighted_hist(values, weights, bins)
        histograms.append((label, hist, err, color))

    global_scale = None
    if normalization == "global":
        global_scale = np.sum([np.sum(hist, dtype=np.float64) for _, hist, _, _ in histograms], dtype=np.float64)

    plotted = False
    for label, hist, err, color in histograms:
        plotted_hist, plotted_err = normalize_hist_and_error(hist, err, normalization, global_scale)
        stairs = ax.stairs(
            plotted_hist,
            bins,
            linewidth=1.8,
            label=label,
            color=color,
        )
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

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(alpha=0.25)
    if plotted:
        ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def collect_category_series(matrix, feature_index, class_names):
    feature_values = matrix[:, feature_index].astype(np.float64)
    weights = matrix[:, WEIGHT_INDEX].astype(np.float64)
    labels = matrix[:, LABEL_INDEX].astype(np.int64)

    series_list = []
    for label_index, category_name in enumerate(class_names):
        mask = labels == label_index
        if not np.any(mask):
            continue
        series_list.append((category_name, feature_values[mask], weights[mask], get_category_color(category_name)))
    return series_list


def collect_channel_series(matrix, feature_index, channel_id_to_name):
    feature_values = matrix[:, feature_index].astype(np.float64)
    weights = matrix[:, WEIGHT_INDEX].astype(np.float64)
    channel_ids = matrix[:, CHANNEL_INDEX].astype(np.int64)

    series_by_channel = {}
    for channel_id, channel_name in channel_id_to_name.items():
        mask = channel_ids == channel_id
        if not np.any(mask):
            continue
        series_by_channel[channel_name] = (feature_values[mask], weights[mask], get_channel_color(channel_name))
    return series_by_channel


def build_prediction_group_labels(prediction_table, dataset_size, group_by):
    group_labels = np.full(dataset_size, -1, dtype=np.int64)
    row_index = prediction_table["rowIndex"].astype(np.int64)
    if np.any(row_index < 0) or np.any(row_index >= dataset_size):
        raise RuntimeError("Prediction rowIndex is out of range for the provided dataset.")
    group_labels[row_index] = prediction_table[group_by].astype(np.int64)
    return group_labels


def plot_focus_channel_by_prediction(matrix, prediction_table, channel_id_to_name, output_dir, nbins, class_names, focus_channel, group_by):
    matching_channel_ids = [
        channel_id for channel_id, channel_name in channel_id_to_name.items() if channel_name == focus_channel
    ]
    if not matching_channel_ids:
        raise RuntimeError(f"Unable to find channel '{focus_channel}' in channel_map.json.")

    group_labels = build_prediction_group_labels(prediction_table, matrix.shape[0], group_by)
    channel_ids = matrix[:, CHANNEL_INDEX].astype(np.int64)
    channel_mask = np.isin(channel_ids, matching_channel_ids)
    if not np.any(channel_mask):
        raise RuntimeError(f"No events from channel '{focus_channel}' are present in the dataset.")

    valid_mask = channel_mask & (group_labels >= 0)
    if not np.any(valid_mask):
        raise RuntimeError(f"No matched prediction rows were found for channel '{focus_channel}'.")

    target_dir = output_dir / "channel_by_prediction" / f"{focus_channel}_{group_by}"
    target_dir.mkdir(parents=True, exist_ok=True)
    weights = matrix[:, WEIGHT_INDEX].astype(np.float64)

    for feature_index, feature_name in enumerate(FEATURE_FIELD_NAMES):
        feature_values = matrix[:, feature_index].astype(np.float64)
        series_list = []
        for class_index, class_name in enumerate(class_names):
            mask = valid_mask & (group_labels == class_index)
            if not np.any(mask):
                continue
            series_list.append(
                (
                    f"{focus_channel} -> {class_name}",
                    feature_values[mask],
                    weights[mask],
                    get_category_color(class_name),
                )
            )
        if not series_list:
            continue
        bins = build_feature_bins(feature_name, nbins)
        plot_hist_series(
            target_dir / f"{feature_name}.png",
            title=f"{feature_name} in {focus_channel} grouped by {group_by}",
            xlabel=feature_name,
            ylabel="Normalized weighted entries",
            bins=bins,
            series_list=series_list,
            normalization="per_hist",
        )


def plot_category_comparison(matrix, output_dir, nbins, class_names):
    target_dir = output_dir / "category_comparison"
    target_dir.mkdir(parents=True, exist_ok=True)

    for feature_index, feature_name in enumerate(FEATURE_FIELD_NAMES):
        series_list = collect_category_series(matrix, feature_index, class_names)
        bins = build_feature_bins(feature_name, nbins)
        plot_hist_series(
            target_dir / f"{feature_name}.png",
            title=f"{feature_name} by category",
            xlabel=feature_name,
            ylabel="Normalized weighted entries",
            bins=bins,
            series_list=series_list,
            normalization="per_hist",
        )


def plot_category_internal(matrix, channel_id_to_name, output_dir, nbins, schema):
    target_dir = output_dir / "category_internal"
    global_dir = target_dir / "global_normalized"
    per_hist_dir = target_dir / "per_channel_normalized"
    global_dir.mkdir(parents=True, exist_ok=True)
    per_hist_dir.mkdir(parents=True, exist_ok=True)

    channel_series = {
        feature_index: collect_channel_series(matrix, feature_index, channel_id_to_name)
        for feature_index in range(FEATURE_DIM)
    }

    channel_to_category = schema.channel_to_class
    channel_order = {channel.channel_name: index for index, channel in enumerate(schema.channels)}

    for feature_index, feature_name in enumerate(FEATURE_FIELD_NAMES):
        bins = build_feature_bins(feature_name, nbins)

        for category_name in schema.class_names:
            series_list = [
                (channel_name, values, weights, color)
                for channel_name, (values, weights, color) in sorted(
                    channel_series[feature_index].items(),
                    key=lambda item: channel_order.get(item[0], 999),
                )
                if channel_to_category.get(channel_name) == category_name
            ]
            if not series_list:
                continue

            plot_hist_series(
                global_dir / f"{category_name}_{feature_name}_global.png",
                title=f"{feature_name} in {category_name}",
                xlabel=feature_name,
                ylabel="Weighted entries / total category weight",
                bins=bins,
                series_list=series_list,
                normalization="global",
            )
            plot_hist_series(
                per_hist_dir / f"{category_name}_{feature_name}_per_channel.png",
                title=f"{feature_name} in {category_name}",
                xlabel=feature_name,
                ylabel="Normalized weighted entries",
                bins=bins,
                series_list=series_list,
                normalization="per_hist",
            )


def main():
    args = parse_args()
    schema = get_signal_schema(args.class_schema)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    matrix = load_dataset(args.dataset)
    channel_id_to_name = load_channel_id_to_name(args.channel_map, schema)

    if not channel_id_to_name:
        raise RuntimeError("No valid channels found in channel_map.json.")

    plot_category_comparison(matrix, output_dir, args.nbins, list(schema.class_names))
    plot_category_internal(matrix, channel_id_to_name, output_dir, args.nbins, schema)
    if args.predictions is not None or args.focus_channel is not None:
        if args.predictions is None or args.focus_channel is None:
            raise RuntimeError("--predictions and --focus-channel must be provided together.")
        prediction_table = load_prediction_table(args.predictions, list(schema.class_names))
        plot_focus_channel_by_prediction(
            matrix,
            prediction_table,
            channel_id_to_name,
            output_dir,
            args.nbins,
            list(schema.class_names),
            args.focus_channel,
            args.group_by,
        )

    total_plots = len(FEATURE_FIELD_NAMES) + len(FEATURE_FIELD_NAMES) * len(schema.class_names) * 2
    print(f"saved output dir: {output_dir}")
    print(f"saved plots: {total_plots}")


if __name__ == "__main__":
    main()
