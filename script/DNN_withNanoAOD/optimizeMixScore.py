#!/usr/bin/env python3
import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_XS_JSON = SCRIPT_DIR.parent.parent / "json" / "XS" / "Run3.json"
DEFAULT_LUMI_JSON = SCRIPT_DIR.parent.parent / "json" / "Lumi" / "Run3.json"
DEFAULT_SIGNAL_PATTERNS = ("TTHH",)
DEFAULT_SCORE_PREFIX = "Signal"
DEFAULT_CLASS_NAMES = ("TTHH", "tt_b", "ttX_like", "TTTT")
DEFAULT_TREE_NAME = "Events"
DEFAULT_WEIGHT_EXPR = "genWeight"
DEFAULT_EPSILON = 1e-6
DEFAULT_SCORE_MODE = "bound"


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Optimize a signal-vs-background mix score using simplex grid scanning. "
            "Signal is taken from TTHH* channels, while every channel in the background json is "
            "treated as background and normalized with XS/genWeightSum."
        )
    )
    parser.add_argument("--signal-json", required=True, help="Sample json containing signal channels.")
    parser.add_argument("--background-json", required=True, help="Sample json containing background channels.")
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory for optimization summaries and tables.",
    )
    parser.add_argument(
        "--xs-json",
        default=str(DEFAULT_XS_JSON),
        help="Cross section json path.",
    )
    parser.add_argument(
        "--lumi-json",
        default=str(DEFAULT_LUMI_JSON),
        help="Integrated luminosity json path.",
    )
    parser.add_argument(
        "--total-lumi",
        type=float,
        default=None,
        help="Optional explicit integrated luminosity in fb^-1. If omitted, infer from lumi json.",
    )
    parser.add_argument(
        "--tree-name",
        default=DEFAULT_TREE_NAME,
        help="Tree name to read from ROOT files.",
    )
    parser.add_argument(
        "--score-prefix",
        default=DEFAULT_SCORE_PREFIX,
        help="Prediction branch prefix written by writeSignalPredictionToRoot.py.",
    )
    parser.add_argument(
        "--score-fields",
        nargs="+",
        default=list(DEFAULT_CLASS_NAMES),
        metavar="FIELD",
        help=(
            "Structured score field suffixes in order: signal, background1, background2, ... "
            "Branch names are built as <prefix>_score_<field>. Supports three-class and four-class outputs."
        ),
    )
    parser.add_argument(
        "--score-mode",
        choices=["bound", "unbound"],
        default=DEFAULT_SCORE_MODE,
        help=(
            "Mix-score definition. "
            "'bound': w_sig*S / (w_sig*S + sum_i w_i*B_i), with all weights summing to 1. "
            "'unbound': S / (sum_i w_i*B_i), with only background weights summing to 1."
        ),
    )
    parser.add_argument(
        "--signal-patterns",
        nargs="+",
        default=list(DEFAULT_SIGNAL_PATTERNS),
        help="Channel-name substrings treated as signal. Defaults to TTHH.",
    )
    parser.add_argument(
        "--weight-branch",
        default=DEFAULT_WEIGHT_EXPR,
        help="Per-event weight branch already stored in the ROOT files.",
    )
    parser.add_argument(
        "--ttbb-scale",
        type=float,
        default=1.0,
        help="Optional extra normalization factor applied only to TTBB* channels.",
    )
    parser.add_argument(
        "--step",
        type=float,
        default=0.1,
        help="Simplex grid step size for the score weights.",
    )
    parser.add_argument(
        "--min-signal-weight",
        type=float,
        default=None,
        help=(
            "Minimum allowed weight for the signal score in the bounded mix. "
            "Defaults to --step, which removes the degenerate w_signal=0 solution."
        ),
    )
    parser.add_argument(
        "--n-cut-scan",
        type=int,
        default=200,
        help="Number of approximately uniform ranked cut positions tested per weight point.",
    )
    parser.add_argument(
        "--metric",
        choices=["asimov", "s_over_sqrt_b"],
        default="asimov",
        help="Counting significance metric used for optimization.",
    )
    parser.add_argument(
        "--epsilon",
        type=float,
        default=DEFAULT_EPSILON,
        help="Small positive number added to the denominator for numerical stability.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print progress information while loading and optimizing.",
    )
    return parser.parse_args()


def log_verbose(enabled, message):
    if enabled:
        print(f"[optimizeMixScore] {message}")


def load_json(path):
    with open(path, "r") as handle:
        return json.load(handle)


def infer_total_lumi(lumi_info, signal_json_path, background_json_path):
    signal_stem = Path(signal_json_path).stem
    background_stem = Path(background_json_path).stem
    stem_text = f"{signal_stem} {background_stem}"
    if "2024" in stem_text:
        target_prefix = "Run2024"
    elif "2023" in stem_text:
        target_prefix = "Run2023"
    elif "2022" in stem_text:
        target_prefix = "Run2022"
    else:
        target_prefix = None

    if target_prefix is not None:
        matched = {key: value for key, value in lumi_info.items() if str(key).startswith(target_prefix)}
        if matched:
            return float(sum(matched.values())), sorted(matched.keys())

    return float(sum(lumi_info.values())), sorted(lumi_info.keys())


def channel_matches_signal(channel_name, signal_patterns):
    channel_upper = str(channel_name).upper()
    return any(pattern.upper() in channel_upper for pattern in signal_patterns)


def channel_extra_scale(channel_name, ttbb_scale):
    if str(channel_name).startswith("TTBB"):
        return float(ttbb_scale)
    return 1.0


def build_channel_file_map(sample_dict):
    channel_files = {}
    for channel, directory in sample_dict["dir"].items():
        file_list = sample_dict["file"].get(channel, [])
        channel_files[channel] = [directory + fname for fname in file_list]
    return channel_files


def read_genweightsum(file_path):
    try:
        with uproot.open(file_path) as root_file:
            if "genWeightSum" not in root_file:
                return None
            obj = root_file["genWeightSum"]
            if hasattr(obj, "values"):
                values = obj.values()
                return float(np.sum(values, dtype=np.float64))
            return None
    except Exception:
        return None


def collect_valid_files_and_sumw(file_list, tree_name, verbose=False, label=""):
    valid_files = []
    total_sumw = 0.0
    for file_path in file_list:
        sumw = read_genweightsum(file_path)
        if sumw is None:
            log_verbose(verbose, f"skip file without usable genWeightSum ({label}): {file_path}")
            continue
        try:
            with uproot.open(file_path) as root_file:
                if tree_name not in root_file and "Events" not in root_file:
                    log_verbose(verbose, f"skip file without Events tree ({label}): {file_path}")
                    continue
        except Exception:
            log_verbose(verbose, f"skip unreadable ROOT file ({label}): {file_path}")
            continue
        valid_files.append(file_path)
        total_sumw += sumw
    return valid_files, total_sumw


def load_channel_arrays(file_list, tree_name, branch_names, norm_factor, verbose=False, label=""):
    arrays = []
    for file_path in file_list:
        try:
            with uproot.open(file_path) as root_file:
                resolved_tree = tree_name if tree_name in root_file else "Events"
                tree = root_file[resolved_tree]
                content = tree.arrays(branch_names, library="np")
        except Exception as exc:
            log_verbose(verbose, f"skip unreadable tree ({label}): {file_path} ({exc})")
            continue

        if not content:
            continue
        event_count = len(content[branch_names[0]])
        if event_count == 0:
            continue

        channel_arrays = {name: np.asarray(content[name], dtype=np.float64) for name in branch_names}
        channel_arrays["scaled_weight"] = channel_arrays[branch_names[-1]] * norm_factor
        arrays.append(channel_arrays)
    return arrays


def concatenate_channel_payload(channel_payloads):
    if not channel_payloads:
        return None
    merged = {}
    for key in channel_payloads[0]:
        merged[key] = np.concatenate([payload[key] for payload in channel_payloads], axis=0)
    return merged


def asimov_significance(signal_yield, background_yield):
    if signal_yield <= 0.0 or background_yield <= 0.0:
        return 0.0
    return math.sqrt(
        2.0 * ((signal_yield + background_yield) * math.log1p(signal_yield / background_yield) - signal_yield)
    )


def simple_significance(signal_yield, background_yield):
    if signal_yield <= 0.0 or background_yield <= 0.0:
        return 0.0
    return signal_yield / math.sqrt(background_yield)


def build_simplex_grid(step, n_weights):
    if step <= 0.0 or step > 1.0:
        raise RuntimeError(f"step must be in (0, 1], got {step}")
    if n_weights < 2:
        raise RuntimeError(f"At least signal and one background score are required, got {n_weights}")
    inv = round(1.0 / step)
    if not np.isclose(inv * step, 1.0, atol=1e-8):
        raise RuntimeError("step must divide 1 exactly, e.g. 0.5, 0.25, 0.2, 0.1, 0.05")

    grid = []

    def fill(prefix, remaining_units, slots_left):
        if slots_left == 1:
            grid.append(tuple(round(value * step, 10) for value in [*prefix, remaining_units]))
            return
        for units in range(remaining_units + 1):
            fill([*prefix, units], remaining_units - units, slots_left - 1)

    fill([], inv, n_weights)
    return grid


def build_scan_indices(n_events, n_cut_scan):
    if n_events <= 0:
        return np.asarray([], dtype=np.int64)
    count = min(max(1, n_cut_scan), n_events)
    indices = np.linspace(1, n_events, num=count, dtype=np.int64) - 1
    return np.unique(indices)


def get_weight_field_names(score_fields, score_mode):
    if score_mode == "bound":
        return [f"w_{field}" for field in score_fields]
    return [f"w_{field}" for field in score_fields[1:]]


def build_candidate_weight_grid(args):
    n_score_fields = len(args.score_fields)
    if args.score_mode == "bound":
        min_signal_weight = args.step if args.min_signal_weight is None else float(args.min_signal_weight)
        if min_signal_weight < 0.0 or min_signal_weight > 1.0:
            raise RuntimeError(f"--min-signal-weight must be in [0, 1], got {min_signal_weight}")
        weight_grid_all = build_simplex_grid(args.step, n_score_fields)
        weight_grid = [weights for weights in weight_grid_all if weights[0] + 1e-12 >= min_signal_weight]
        if not weight_grid:
            raise RuntimeError(
                f"No simplex weight points left after requiring w_signal >= {min_signal_weight}. "
                "Use a smaller --min-signal-weight or a coarser --step."
            )
        return weight_grid_all, weight_grid, min_signal_weight

    if args.min_signal_weight is not None:
        raise RuntimeError("--min-signal-weight is only valid in --score-mode bound.")
    background_weight_grid = build_simplex_grid(args.step, n_score_fields - 1)
    return background_weight_grid, background_weight_grid, None


def evaluate_ranked_score(score_values, signal_weights, background_weights, candidate_indices, metric_fn):
    order = np.argsort(score_values)[::-1]
    ordered_scores = score_values[order]
    ordered_signal = signal_weights[order]
    ordered_background = background_weights[order]
    signal_cumsum = np.cumsum(ordered_signal, dtype=np.float64)
    background_cumsum = np.cumsum(ordered_background, dtype=np.float64)

    best = None
    for idx in candidate_indices:
        signal_yield = float(signal_cumsum[idx])
        background_yield = float(background_cumsum[idx])
        if signal_yield <= 0.0 or background_yield <= 0.0:
            continue
        significance = float(metric_fn(signal_yield, background_yield))
        threshold = float(ordered_scores[idx])
        record = {
            "cut_index": int(idx),
            "threshold": threshold,
            "signal_yield": signal_yield,
            "background_yield": background_yield,
            "significance": significance,
        }
        if best is None or significance > best["significance"]:
            best = record
    return best


def build_ranked_curve(score_values, signal_weights, background_weights, candidate_indices, metric_fn):
    order = np.argsort(score_values)[::-1]
    ordered_scores = score_values[order]
    ordered_signal = signal_weights[order]
    ordered_background = background_weights[order]
    signal_cumsum = np.cumsum(ordered_signal, dtype=np.float64)
    background_cumsum = np.cumsum(ordered_background, dtype=np.float64)

    curve = []
    for idx in candidate_indices:
        signal_yield = float(signal_cumsum[idx])
        background_yield = float(background_cumsum[idx])
        significance = float(metric_fn(signal_yield, background_yield))
        curve.append(
            {
                "cut_index": int(idx),
                "threshold": float(ordered_scores[idx]),
                "signal_yield": signal_yield,
                "background_yield": background_yield,
                "significance": significance,
            }
        )
    return curve


def compute_mix_values(scores, weights, score_mode, epsilon):
    scores = np.asarray(scores, dtype=np.float64)
    weights_array = np.asarray(weights, dtype=np.float64)
    signal_score = scores[:, 0]
    background_scores = scores[:, 1:]

    if background_scores.shape[1] != len(weights_array) - (1 if score_mode == "bound" else 0):
        raise RuntimeError(
            f"Weight dimension does not match score mode '{score_mode}': "
            f"scores have {scores.shape[1]} classes but weights are {weights_array}."
        )

    if score_mode == "bound":
        denominator = np.dot(scores, weights_array) + epsilon
        numerator = weights_array[0] * signal_score
        return numerator / denominator

    denominator = np.dot(background_scores, weights_array) + epsilon
    return signal_score / denominator


def optimize_mix(signal_payload, background_payload, candidate_weights, candidate_indices, metric_fn, score_mode, epsilon, verbose=False):
    signal_scores = signal_payload["scores"]
    signal_weights = signal_payload["weight"]

    background_scores = background_payload["scores"]
    background_weights = background_payload["weight"]

    ranking = []
    for idx, weights in enumerate(candidate_weights, start=1):
        signal_mix = compute_mix_values(signal_scores, weights, score_mode, epsilon)
        background_mix = compute_mix_values(background_scores, weights, score_mode, epsilon)

        all_scores = np.concatenate([signal_mix, background_mix], axis=0)
        all_signal_weights = np.concatenate([signal_weights, np.zeros_like(background_weights)], axis=0)
        all_background_weights = np.concatenate([np.zeros_like(signal_weights), background_weights], axis=0)

        best = evaluate_ranked_score(all_scores, all_signal_weights, all_background_weights, candidate_indices, metric_fn)
        if best is None:
            continue
        best["weights"] = [float(weight) for weight in weights]
        ranking.append(best)
        if verbose and (idx == 1 or idx % 10 == 0 or idx == len(candidate_weights)):
            log_verbose(
                True,
                f"tested {idx}/{len(candidate_weights)} weight points; current best Z={max(r['significance'] for r in ranking):.6g}",
            )

    ranking.sort(key=lambda item: item["significance"], reverse=True)
    return ranking


def build_mix_scores(signal_payload, background_payload, weights, score_mode, epsilon):
    signal_mix = compute_mix_values(signal_payload["scores"], weights, score_mode, epsilon)
    background_mix = compute_mix_values(background_payload["scores"], weights, score_mode, epsilon)
    return np.concatenate([signal_mix, background_mix], axis=0)


def write_curve_csv(path, curve_map):
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for label, curve in curve_map.items():
        for point in curve:
            rows.append(
                {
                    "label": label,
                    "cut_index": point["cut_index"],
                    "threshold": point["threshold"],
                    "signal_yield": point["signal_yield"],
                    "background_yield": point["background_yield"],
                    "significance": point["significance"],
                }
            )
    dump_csv(path, rows)


def plot_significance_curves(path, curve_map, metric_name):
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(9, 6))
    for label, curve in curve_map.items():
        thresholds = [point["threshold"] for point in curve]
        significances = [point["significance"] for point in curve]
        plt.plot(thresholds, significances, linewidth=2, label=label)
    plt.xlabel("Score Cut")
    plt.ylabel(f"{metric_name} significance")
    plt.title("Significance vs cut")
    plt.legend(fontsize=9)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(path, dpi=160)
    plt.close()


def plot_best_significance_points(path, point_rows, metric_name):
    path.parent.mkdir(parents=True, exist_ok=True)
    labels = [row["label"] for row in point_rows]
    values = [row["significance"] for row in point_rows]
    x = np.arange(len(labels), dtype=np.float64)

    plt.figure(figsize=(10, 6))
    plt.scatter(x, values, s=80)
    plt.plot(x, values, linewidth=1.5, alpha=0.7)
    plt.xticks(x, labels, rotation=20, ha="right")
    plt.ylabel(f"best {metric_name} significance")
    plt.xlabel("score combination")
    plt.title("Best significance by score combination")
    plt.grid(True, axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(path, dpi=160)
    plt.close()


def dump_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as handle:
        json.dump(payload, handle, indent=2)


def dump_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def summarize_channel_yields(channel_summaries):
    out = []
    for item in channel_summaries:
        out.append(
            {
                "channel": item["channel"],
                "files": item["n_files"],
                "genWeightSum": item["genweightsum"],
                "norm_factor": item["norm_factor"],
                "event_count": item["event_count"],
                "weighted_yield": item["weighted_yield"],
            }
        )
    return out


def format_weights(weights, precision=6):
    return "(" + ", ".join(f"{weight:.{precision}g}" for weight in weights) + ")"


def format_weights_fixed(weights, precision=2):
    return "(" + ",".join(f"{weight:.{precision}f}" for weight in weights) + ")"


def load_labeled_samples(
    sample_json_path,
    xs_info,
    tree_name,
    branch_names,
    total_lumi,
    signal_patterns,
    take_signal,
    ttbb_scale=1.0,
    verbose=False,
):
    sample_dict = load_json(sample_json_path)
    channel_files = build_channel_file_map(sample_dict)
    payloads = []
    summaries = []

    for channel, file_list in channel_files.items():
        is_signal_channel = channel_matches_signal(channel, signal_patterns)
        if take_signal != is_signal_channel:
            continue
        if channel not in xs_info:
            log_verbose(verbose, f"skip channel without XS: {channel}")
            continue
        valid_files, total_sumw = collect_valid_files_and_sumw(file_list, tree_name, verbose, channel)
        if not valid_files or total_sumw <= 0.0:
            log_verbose(verbose, f"skip channel with no valid files or non-positive genWeightSum: {channel}")
            continue
        extra_scale = channel_extra_scale(channel, ttbb_scale)
        norm_factor = float(xs_info[channel]) * float(total_lumi) * 1000.0 * extra_scale / float(total_sumw)
        channel_arrays = load_channel_arrays(valid_files, tree_name, branch_names, norm_factor, verbose, channel)
        merged = concatenate_channel_payload(channel_arrays)
        if merged is None:
            continue
        score_matrix = np.column_stack([merged[name] for name in branch_names[:-1]])
        payloads.append(
            {
                "scores": score_matrix,
                "sig": score_matrix[:, 0],
                "weight": merged["scaled_weight"],
            }
        )
        summaries.append(
            {
                "channel": channel,
                "n_files": len(valid_files),
                "genweightsum": float(total_sumw),
                "extra_scale": float(extra_scale),
                "norm_factor": float(norm_factor),
                "event_count": int(len(merged["scaled_weight"])),
                "weighted_yield": float(np.sum(merged["scaled_weight"], dtype=np.float64)),
            }
        )
        log_verbose(
            verbose,
            f"loaded channel {channel}: files={len(valid_files)}, events={len(merged['scaled_weight'])}, "
            f"extraScale={extra_scale:.6g}, yield={np.sum(merged['scaled_weight'], dtype=np.float64):.6g}",
        )

    merged = concatenate_channel_payload(payloads)
    if merged is None:
        return None, summaries
    return merged, summaries


def main():
    args = parse_args()
    if len(args.score_fields) < 2:
        raise RuntimeError("--score-fields must contain at least signal and one background score field.")
    if args.score_mode == "bound" and len(args.score_fields) < 3:
        raise RuntimeError("--score-mode bound expects at least one signal score and two total classes.")
    if args.score_mode == "unbound" and len(args.score_fields) < 3:
        raise RuntimeError("--score-mode unbound expects one signal score and at least two background score fields.")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    xs_info = load_json(args.xs_json)
    lumi_info = load_json(args.lumi_json)
    if args.total_lumi is None:
        total_lumi, lumi_keys = infer_total_lumi(lumi_info, args.signal_json, args.background_json)
    else:
        total_lumi = float(args.total_lumi)
        lumi_keys = ["explicit"]

    score_branches = [f"{args.score_prefix}_score_{field}" for field in args.score_fields]
    branch_names = score_branches + [args.weight_branch]
    metric_fn = asimov_significance if args.metric == "asimov" else simple_significance

    log_verbose(args.verbose, f"score branches: {score_branches}")
    log_verbose(args.verbose, f"weight branch: {args.weight_branch}")
    log_verbose(args.verbose, f"using total lumi = {total_lumi:.6g} fb^-1 from {lumi_keys}")

    signal_payload, signal_channel_summaries = load_labeled_samples(
        args.signal_json,
        xs_info,
        args.tree_name,
        branch_names,
        total_lumi,
        args.signal_patterns,
        take_signal=True,
        ttbb_scale=args.ttbb_scale,
        verbose=args.verbose,
    )
    background_payload, background_channel_summaries = load_labeled_samples(
        args.background_json,
        xs_info,
        args.tree_name,
        branch_names,
        total_lumi,
        args.signal_patterns,
        take_signal=False,
        ttbb_scale=args.ttbb_scale,
        verbose=args.verbose,
    )

    if signal_payload is None:
        raise RuntimeError("No valid signal events were loaded.")
    if background_payload is None:
        raise RuntimeError("No valid background events were loaded.")

    n_signal = len(signal_payload["weight"])
    n_background = len(background_payload["weight"])
    log_verbose(
        args.verbose,
        f"merged payload sizes: signal events={n_signal}, background events={n_background}",
    )

    total_events = n_signal + n_background
    candidate_indices = build_scan_indices(total_events, args.n_cut_scan)
    weight_grid_all, weight_grid, min_signal_weight = build_candidate_weight_grid(args)
    if args.score_mode == "bound":
        log_verbose(
            args.verbose,
            f"simplex grid points: {len(weight_grid)} after w_signal >= {min_signal_weight} "
            f"(from {len(weight_grid_all)} total) with step={args.step}",
        )
    else:
        log_verbose(
            args.verbose,
            f"background-only simplex grid points: {len(weight_grid)} with step={args.step}",
        )
    log_verbose(args.verbose, f"candidate cut count per score: {len(candidate_indices)}")

    baseline_scores = np.concatenate([signal_payload["sig"], background_payload["sig"]], axis=0)
    baseline_signal_weights = np.concatenate(
        [signal_payload["weight"], np.zeros_like(background_payload["weight"])],
        axis=0,
    )
    baseline_background_weights = np.concatenate(
        [np.zeros_like(signal_payload["weight"]), background_payload["weight"]],
        axis=0,
    )
    baseline_best = evaluate_ranked_score(
        baseline_scores,
        baseline_signal_weights,
        baseline_background_weights,
        candidate_indices,
        metric_fn,
    )
    if baseline_best is None:
        raise RuntimeError("Baseline pure signal-score scan found no valid cut with positive S and B.")
    baseline_best["weights"] = None
    baseline_curve = build_ranked_curve(
        baseline_scores,
        baseline_signal_weights,
        baseline_background_weights,
        candidate_indices,
        metric_fn,
    )

    ranking = optimize_mix(
        signal_payload,
        background_payload,
        weight_grid,
        candidate_indices,
        metric_fn,
        args.score_mode,
        args.epsilon,
        verbose=args.verbose,
    )
    if not ranking:
        raise RuntimeError("No valid weight point produced a positive S and B significance result.")

    best = ranking[0]
    n_score_fields = len(args.score_fields)
    n_weight_fields = n_score_fields if args.score_mode == "bound" else n_score_fields - 1
    uniform_weights = tuple([1.0 / n_weight_fields] * n_weight_fields)
    uniform_scores = build_mix_scores(signal_payload, background_payload, uniform_weights, args.score_mode, args.epsilon)
    uniform_best = evaluate_ranked_score(
        uniform_scores,
        baseline_signal_weights,
        baseline_background_weights,
        candidate_indices,
        metric_fn,
    )
    uniform_curve = build_ranked_curve(
        uniform_scores,
        baseline_signal_weights,
        baseline_background_weights,
        candidate_indices,
        metric_fn,
    )
    comparison = {
        "baseline_signal_score": baseline_best,
        "uniform_mix_score": uniform_best,
        "best_mix_score": best,
        "relative_gain": (
            best["significance"] / baseline_best["significance"]
            if baseline_best["significance"] > 0.0
            else None
        ),
    }

    summary = {
        "config": {
            "signal_json": args.signal_json,
            "background_json": args.background_json,
            "xs_json": args.xs_json,
            "lumi_json": args.lumi_json,
            "total_lumi": total_lumi,
            "lumi_keys": lumi_keys,
            "tree_name": args.tree_name,
            "score_prefix": args.score_prefix,
            "score_fields": args.score_fields,
            "score_mode": args.score_mode,
            "score_branches": score_branches,
            "weight_branch": args.weight_branch,
            "ttbb_scale": args.ttbb_scale,
            "step": args.step,
            "min_signal_weight": min_signal_weight,
            "n_cut_scan": args.n_cut_scan,
            "metric": args.metric,
            "epsilon": args.epsilon,
            "signal_patterns": args.signal_patterns,
        },
        "event_summary": {
            "signal_events": int(n_signal),
            "background_events": int(n_background),
            "signal_weighted_yield_total": float(np.sum(signal_payload["weight"], dtype=np.float64)),
            "background_weighted_yield_total": float(np.sum(background_payload["weight"], dtype=np.float64)),
        },
        "signal_channels": summarize_channel_yields(signal_channel_summaries),
        "background_channels": summarize_channel_yields(background_channel_summaries),
        "comparison": comparison,
        "top_weight_points": ranking[:20],
    }

    csv_rows = []
    weight_field_names = get_weight_field_names(args.score_fields, args.score_mode)
    for item in ranking:
        row = {weight_field_names[idx]: item["weights"][idx] for idx in range(len(weight_field_names))}
        row.update(
            {
                "threshold": item["threshold"],
                "cut_index": item["cut_index"],
                "signal_yield": item["signal_yield"],
                "background_yield": item["background_yield"],
                "significance": item["significance"],
            }
        )
        csv_rows.append(row)

    curve_map = {
        "pure_signal": baseline_curve,
        f"uniform_{args.score_mode}_mix_1over{n_weight_fields}": uniform_curve,
    }
    best_point_rows = [
        {
            "label": "pure_signal",
            "significance": baseline_best["significance"],
            "threshold": baseline_best["threshold"],
            "signal_yield": baseline_best["signal_yield"],
            "background_yield": baseline_best["background_yield"],
        }
    ]
    if uniform_best is not None:
        best_point_rows.append(
            {
                "label": f"uniform_{args.score_mode}_mix_1over{n_weight_fields}",
                "significance": uniform_best["significance"],
                "threshold": uniform_best["threshold"],
                "signal_yield": uniform_best["signal_yield"],
                "background_yield": uniform_best["background_yield"],
            }
        )
    top_curve_limit = min(5, len(ranking))
    for index in range(top_curve_limit):
        item = ranking[index]
        label = (
            f"top{index + 1}_w="
            f"{format_weights_fixed(item['weights'])}"
        )
        score_values = build_mix_scores(signal_payload, background_payload, item["weights"], args.score_mode, args.epsilon)
        curve_map[label] = build_ranked_curve(
            score_values,
            baseline_signal_weights,
            baseline_background_weights,
            candidate_indices,
            metric_fn,
        )
        best_point_rows.append(
            {
                "label": label,
                "significance": item["significance"],
                "threshold": item["threshold"],
                "signal_yield": item["signal_yield"],
                "background_yield": item["background_yield"],
            }
        )

    dump_json(output_dir / "optimization_summary.json", summary)
    dump_csv(output_dir / "weight_scan_ranking.csv", csv_rows)
    write_curve_csv(output_dir / "significance_curves.csv", curve_map)
    dump_csv(output_dir / "best_significance_points.csv", best_point_rows)
    plot_significance_curves(output_dir / "significance_curves.png", curve_map, args.metric)
    plot_best_significance_points(output_dir / "best_significance_points.png", best_point_rows, args.metric)

    print("Optimization finished")
    print(f"metric: {args.metric}")
    print(f"score mode: {args.score_mode}")
    print(f"total lumi [fb^-1]: {total_lumi:.6g}")
    print(
        "baseline signal-score best: "
        f"Z={baseline_best['significance']:.6g}, "
        f"cut={baseline_best['threshold']:.6g}, "
        f"S={baseline_best['signal_yield']:.6g}, "
        f"B={baseline_best['background_yield']:.6g}"
    )
    print(
        "best mix-score point: "
        f"w={format_weights(best['weights'])}, "
        f"Z={best['significance']:.6g}, "
        f"cut={best['threshold']:.6g}, "
        f"S={best['signal_yield']:.6g}, "
        f"B={best['background_yield']:.6g}"
    )
    if comparison["relative_gain"] is not None:
        print(f"relative gain over pure signal score: {comparison['relative_gain']:.6g}")
    print(f"summary json: {output_dir / 'optimization_summary.json'}")
    print(f"scan csv: {output_dir / 'weight_scan_ranking.csv'}")
    print(f"curve csv: {output_dir / 'significance_curves.csv'}")
    print(f"curve plot: {output_dir / 'significance_curves.png'}")
    print(f"best-point csv: {output_dir / 'best_significance_points.csv'}")
    print(f"best-point plot: {output_dir / 'best_significance_points.png'}")


if __name__ == "__main__":
    main()
