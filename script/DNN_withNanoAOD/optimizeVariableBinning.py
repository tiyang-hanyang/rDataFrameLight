import argparse
import csv
import fnmatch
import json
import keyword
import re
from pathlib import Path, PurePosixPath

import numpy as np
import uproot


def parse_args():
    source_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(
        description=(
            "Build variable-width bin edges that flatten a chosen background distribution "
            "after summing weighted channels from a sample json."
        )
    )
    parser.add_argument("--sample-json", "--sampleJson", dest="sample_json", required=True)
    parser.add_argument(
        "--channels",
        nargs="+",
        required=True,
        help="One or more channel names or fnmatch patterns to combine.",
    )
    parser.add_argument("--variable", required=True, help="Branch name or output variable name to optimize.")
    parser.add_argument(
        "--variable-expr",
        default=None,
        help="Optional expression evaluated from branches, e.g. 'A/(0.2*B+0.8*C+1e-6)'.",
    )
    parser.add_argument(
        "--define-config",
        default=None,
        help="Optional json file containing [action, name, expression] definitions.",
    )
    parser.add_argument("--n-bins", type=int, required=True, help="Requested number of bins.")
    parser.add_argument("--tree-name", default="Events")
    parser.add_argument("--config-key", default=None, help="Key name used in the output varConfig snippet.")
    parser.add_argument("--label", default=None, help="Axis label used in the output varConfig snippet.")
    parser.add_argument(
        "--xs-config",
        default=str(source_root / "json" / "XS" / "Run3.json"),
        help="Path to XS json.",
    )
    parser.add_argument(
        "--lumi-config",
        default=str(source_root / "json" / "Lumi" / "Run3.json"),
        help="Path to Lumi json.",
    )
    parser.add_argument("--lumi-key", default=None, help="Key inside lumi json, e.g. RunIII2024Summer24NanoAODv15.")
    parser.add_argument("--lumi", type=float, default=None, help="Direct lumi value. Overrides --lumi-key if set.")
    parser.add_argument(
        "--weight-mode",
        choices=("user_formula", "standard_mc"),
        default="user_formula",
        help=(
            "Event-weight convention. "
            "user_formula: (genWeight * PUWeight * btag_weight) / (genWeightSum * XS * lumi * 1000). "
            "standard_mc: (genWeight * PUWeight * btag_weight) * XS * lumi * 1000 / genWeightSum."
        ),
    )
    parser.add_argument("--genweight-branch", default="genWeight")
    parser.add_argument("--puweight-branch", default="PUWeight")
    parser.add_argument("--btag-weight-branch", default="btag_weight")
    parser.add_argument(
        "--binning-weight-mode",
        choices=("signed", "abs", "positive"),
        default="abs",
        help=(
            "Weights used only for the equal-content optimization. "
            "The output bin summaries always report signed histogram contents."
        ),
    )
    parser.add_argument("--x-min", type=float, default=None, help="Optional lower edge. Underflow is merged into the first bin.")
    parser.add_argument("--x-max", type=float, default=None, help="Optional upper edge. Overflow is merged into the last bin.")
    parser.add_argument("--output-dir", default="optimize_variable_binning")
    parser.add_argument("--summary-name", default="optimization_summary.json")
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def log_verbose(enabled, message):
    if enabled:
        print(f"[optimizeVariableBinning] {message}")


def load_json_with_comments(path):
    with open(path, "r", encoding="utf-8") as handle:
        content = handle.read()
    try:
        return json.loads(content)
    except json.JSONDecodeError:
        sanitized = re.sub(r"//.*?$", "", content, flags=re.MULTILINE)
        sanitized = re.sub(r"#.*?$", "", sanitized, flags=re.MULTILINE)
        sanitized = re.sub(r",(\s*[}\]])", r"\1", sanitized)
        return json.loads(sanitized)


def join_sample_path(base_dir, file_name):
    base_dir = str(base_dir)
    if "/" in base_dir:
        return str(PurePosixPath(base_dir) / file_name)
    return str(Path(base_dir) / file_name)


def extract_identifiers(expression):
    tokens = re.findall(r"\b[A-Za-z_][A-Za-z0-9_]*\b", expression)
    excluded = {
        "and",
        "or",
        "not",
        "True",
        "False",
        "np",
        "numpy",
        "where",
        "abs",
        "sqrt",
        "log",
        "log10",
        "exp",
        "minimum",
        "maximum",
        "clip",
        "pow",
        "sin",
        "cos",
        "tan",
        "arcsin",
        "arccos",
        "arctan",
    }
    result = []
    for token in tokens:
        if token in excluded or keyword.iskeyword(token):
            continue
        result.append(token)
    return list(dict.fromkeys(result))


def load_define_expressions(path):
    payload = load_json_with_comments(path)
    if not isinstance(payload, list):
        raise RuntimeError(f"Define config must be a json array: {path}")
    definitions = {}
    for item in payload:
        if not isinstance(item, list) or len(item) < 3:
            continue
        action, name, expression = item[0], item[1], item[2]
        if action not in ("define", "redefine"):
            continue
        if not isinstance(name, str) or not isinstance(expression, str):
            continue
        definitions[name] = expression
    return definitions


def resolve_variable_expression(args):
    if args.variable_expr is not None:
        return args.variable_expr, {}
    if args.define_config is None:
        return None, {}
    definitions = load_define_expressions(args.define_config)
    if args.variable not in definitions:
        raise RuntimeError(
            f"Variable {args.variable} was not found in define config {args.define_config}"
        )
    return definitions[args.variable], definitions


def collect_required_branches(variable_name, variable_expression, definitions, weight_branches):
    required = list(weight_branches)
    pending = []
    if variable_expression is None:
        required.append(variable_name)
    else:
        pending.extend(extract_identifiers(variable_expression))

    seen_defined = set()
    while pending:
        token = pending.pop()
        if token in definitions:
            if token in seen_defined:
                continue
            seen_defined.add(token)
            pending.extend(extract_identifiers(definitions[token]))
        else:
            required.append(token)
    return list(dict.fromkeys(required))


def evaluate_expression(expression, branch_arrays, definitions, cache):
    if expression in cache:
        return cache[expression]

    context = {
        "np": np,
        "numpy": np,
        "where": np.where,
        "abs": np.abs,
        "sqrt": np.sqrt,
        "log": np.log,
        "log10": np.log10,
        "exp": np.exp,
        "minimum": np.minimum,
        "maximum": np.maximum,
        "clip": np.clip,
        "pow": np.power,
        "sin": np.sin,
        "cos": np.cos,
        "tan": np.tan,
        "arcsin": np.arcsin,
        "arccos": np.arccos,
        "arctan": np.arctan,
    }
    for name, values in branch_arrays.items():
        context[name] = values

    for name in extract_identifiers(expression):
        if name in definitions and name not in context:
            context[name] = evaluate_expression(definitions[name], branch_arrays, definitions, cache)

    try:
        result = eval(expression, {"__builtins__": {}}, context)
    except NameError as exc:
        raise RuntimeError(f"Failed to resolve expression '{expression}': {exc}") from exc
    except Exception as exc:
        raise RuntimeError(f"Failed to evaluate expression '{expression}': {exc}") from exc

    cache[expression] = np.asarray(result, dtype=np.float64)
    return cache[expression]


def resolve_channel_files(sample_json_path, channel_patterns):
    payload = load_json_with_comments(sample_json_path)
    dir_map = payload.get("dir")
    file_map = payload.get("file")
    if not isinstance(dir_map, dict) or not isinstance(file_map, dict):
        raise RuntimeError(f"Sample json must contain object fields 'dir' and 'file': {sample_json_path}")

    available_channels = sorted(file_map.keys())
    selected_channels = []
    for pattern in channel_patterns:
        matched = [name for name in available_channels if fnmatch.fnmatch(name, pattern)]
        if not matched:
            raise RuntimeError(f"Channel pattern {pattern} did not match anything in {sample_json_path}")
        selected_channels.extend(matched)
    selected_channels = list(dict.fromkeys(selected_channels))

    resolved = []
    for channel in selected_channels:
        if channel not in dir_map or channel not in file_map:
            raise RuntimeError(f"Channel {channel} not found in sample json: {sample_json_path}")
        base_dir = dir_map[channel]
        files = file_map[channel]
        if not files:
            raise RuntimeError(f"Channel {channel} has no files in sample json: {sample_json_path}")
        for one_file in files:
            resolved.append({"channel": channel, "path": join_sample_path(base_dir, one_file)})
    return resolved


def read_genweight_sum(root_file, object_name="genWeightSum"):
    if object_name not in root_file:
        raise RuntimeError(f"Missing {object_name} in ROOT file: {root_file.file.path}")
    hist = root_file[object_name]
    if hasattr(hist, "values"):
        values = np.asarray(hist.values(flow=True), dtype=np.float64)
        return float(np.sum(values, dtype=np.float64))
    if hasattr(hist, "to_numpy"):
        payload = hist.to_numpy(flow=True)
        if isinstance(payload, tuple) and payload:
            values = np.asarray(payload[0], dtype=np.float64)
            return float(np.sum(values, dtype=np.float64))
    raise RuntimeError(f"Unable to read {object_name} from ROOT file: {root_file.file.path}")


def read_lumi(args):
    if args.lumi is not None:
        return float(args.lumi)
    if not args.lumi_key:
        raise RuntimeError("Either --lumi or --lumi-key is required.")
    lumi_map = load_json_with_comments(args.lumi_config)
    if args.lumi_key not in lumi_map:
        raise RuntimeError(f"Lumi key {args.lumi_key} not found in {args.lumi_config}")
    return float(lumi_map[args.lumi_key])


def compute_event_weights(arrays, channel, xs_map, lumi_value, genweight_sum, args):
    if channel not in xs_map:
        raise RuntimeError(f"Channel {channel} not found in XS config {args.xs_config}")
    if genweight_sum == 0.0:
        raise RuntimeError(f"genWeightSum is zero for channel {channel}")

    gen_weight = np.asarray(arrays[args.genweight_branch], dtype=np.float64)
    pu_weight = np.asarray(arrays[args.puweight_branch], dtype=np.float64)
    btag_weight = np.asarray(arrays[args.btag_weight_branch], dtype=np.float64)
    base = gen_weight * pu_weight * btag_weight

    xs_value = float(xs_map[channel])
    if args.weight_mode == "user_formula":
        norm = genweight_sum * xs_value * lumi_value * 1000.0
        if norm == 0.0:
            raise RuntimeError(f"Normalization is zero for channel {channel}")
        return base / norm
    scale = xs_value * lumi_value * 1000.0 / genweight_sum
    return base * scale


def transform_optimization_weights(weights, mode):
    if mode == "signed":
        return np.asarray(weights, dtype=np.float64)
    if mode == "abs":
        return np.abs(np.asarray(weights, dtype=np.float64))
    return np.clip(np.asarray(weights, dtype=np.float64), a_min=0.0, a_max=None)


def clip_for_bounds(values, x_min, x_max):
    work = np.asarray(values, dtype=np.float64)
    if x_min is not None:
        work = np.maximum(work, x_min)
    if x_max is not None:
        work = np.minimum(work, x_max)
    return work


def build_equal_content_edges(values, weights_for_optimization, n_bins, x_min=None, x_max=None):
    if n_bins < 1:
        raise RuntimeError("--n-bins must be positive.")

    clipped_values = clip_for_bounds(values, x_min, x_max)
    finite_mask = np.isfinite(clipped_values) & np.isfinite(weights_for_optimization)
    clipped_values = clipped_values[finite_mask]
    weights_for_optimization = weights_for_optimization[finite_mask]
    if clipped_values.size == 0:
        raise RuntimeError("No finite events remain after filtering.")

    if x_min is None:
        x_min = float(np.min(clipped_values))
    if x_max is None:
        x_max = float(np.max(clipped_values))
    if not x_max > x_min:
        raise RuntimeError(f"Invalid range: x_min={x_min}, x_max={x_max}")

    positive_total = float(np.sum(weights_for_optimization, dtype=np.float64))
    if positive_total <= 0.0:
        raise RuntimeError(
            "Total optimization weight is not positive. "
            "Try --binning-weight-mode abs or positive if the signed sum is unstable."
        )

    order = np.argsort(clipped_values, kind="mergesort")
    sorted_values = clipped_values[order]
    sorted_weights = weights_for_optimization[order]
    if np.any(sorted_weights < 0.0):
        raise RuntimeError(
            "Optimization weights contain negative values after the selected transformation. "
            "Use --binning-weight-mode abs or positive for stable equal-content binning."
        )

    unique_values, inverse = np.unique(sorted_values, return_inverse=True)
    unique_weights = np.zeros(unique_values.shape[0], dtype=np.float64)
    np.add.at(unique_weights, inverse, sorted_weights)

    if unique_values.size < n_bins:
        raise RuntimeError(
            f"Only {unique_values.size} distinct values remain after clipping, "
            f"which is fewer than the requested {n_bins} bins."
        )

    cumulative = np.cumsum(unique_weights, dtype=np.float64)
    targets = positive_total * (np.arange(1, n_bins, dtype=np.float64) / n_bins)

    edges = [float(x_min)]
    last_edge = float(x_min)
    for target in targets:
        idx = int(np.searchsorted(cumulative, target, side="left"))
        idx = min(max(idx, 1), unique_values.size - 1)
        edge = 0.5 * (float(unique_values[idx - 1]) + float(unique_values[idx]))
        if edge <= last_edge:
            edge = float(np.nextafter(last_edge, float(x_max)))
        if edge >= x_max:
            raise RuntimeError(
                "Failed to construct strictly increasing variable bin edges within the requested range. "
                "Try fewer bins or a wider x-range."
            )
        edges.append(edge)
        last_edge = edge
    edges.append(float(x_max))
    return np.asarray(edges, dtype=np.float64), clipped_values, finite_mask


def histogram_summary(values, signed_weights, optimization_weights, edges):
    signed_contents, _ = np.histogram(values, bins=edges, weights=signed_weights)
    abs_contents, _ = np.histogram(values, bins=edges, weights=np.abs(signed_weights))
    opt_contents, _ = np.histogram(values, bins=edges, weights=optimization_weights)
    counts, _ = np.histogram(values, bins=edges)
    rows = []
    for idx in range(len(edges) - 1):
        rows.append(
            {
                "bin": idx + 1,
                "low": float(edges[idx]),
                "high": float(edges[idx + 1]),
                "entries": int(counts[idx]),
                "signed_weight_sum": float(signed_contents[idx]),
                "abs_weight_sum": float(abs_contents[idx]),
                "optimization_weight_sum": float(opt_contents[idx]),
            }
        )
    return rows


def make_variable_config(config_key, label, n_bins, edges):
    return {
        config_key: {
            "label": label,
            "nBins": str(int(n_bins)),
            "binning": [float(f"{edge:.12g}") for edge in edges],
            "binLabels": [],
        }
    }


def save_json(path, payload):
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=4)
        handle.write("\n")


def save_csv(path, rows):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "bin",
                "low",
                "high",
                "entries",
                "signed_weight_sum",
                "abs_weight_sum",
                "optimization_weight_sum",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()
    if (args.x_min is None) != (args.x_max is None):
        raise RuntimeError("--x-min and --x-max must be provided together.")
    if args.n_bins < 1:
        raise RuntimeError("--n-bins must be positive.")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    config_key = args.config_key or args.variable
    label = args.label or args.variable
    variable_expression, definitions = resolve_variable_expression(args)

    xs_map = load_json_with_comments(args.xs_config)
    lumi_value = read_lumi(args)
    resolved_files = resolve_channel_files(args.sample_json, args.channels)

    weight_branches = [
        args.genweight_branch,
        args.puweight_branch,
        args.btag_weight_branch,
    ]
    branch_names = collect_required_branches(args.variable, variable_expression, definitions, weight_branches)

    merged_values = []
    merged_signed_weights = []
    channel_weight_summaries = {}

    for item in resolved_files:
        channel = item["channel"]
        file_path = item["path"]
        log_verbose(args.verbose, f"reading channel={channel} file={file_path}")
        with uproot.open(file_path) as root_file:
            if args.tree_name not in root_file:
                raise RuntimeError(f"Tree {args.tree_name} not found in file: {file_path}")
            tree = root_file[args.tree_name]
            missing_branches = [name for name in branch_names if name not in tree.keys()]
            if missing_branches:
                raise RuntimeError(
                    f"Missing branches in {file_path}: {', '.join(missing_branches)}"
                )
            arrays = tree.arrays(branch_names, library="np")
            genweight_sum = read_genweight_sum(root_file)
            signed_weights = compute_event_weights(arrays, channel, xs_map, lumi_value, genweight_sum, args)
            if variable_expression is None:
                values = np.asarray(arrays[args.variable], dtype=np.float64)
            else:
                values = evaluate_expression(variable_expression, arrays, definitions, {})

        valid_mask = np.isfinite(values) & np.isfinite(signed_weights)
        values = values[valid_mask]
        signed_weights = signed_weights[valid_mask]
        if values.size == 0:
            log_verbose(args.verbose, f"skipping empty finite payload for channel={channel} file={file_path}")
            continue

        merged_values.append(values)
        merged_signed_weights.append(signed_weights)
        channel_summary = channel_weight_summaries.setdefault(
            channel,
            {"entries": 0, "signed_weight_sum": 0.0, "abs_weight_sum": 0.0},
        )
        channel_summary["entries"] += int(values.size)
        channel_summary["signed_weight_sum"] += float(np.sum(signed_weights, dtype=np.float64))
        channel_summary["abs_weight_sum"] += float(np.sum(np.abs(signed_weights), dtype=np.float64))

    if not merged_values:
        raise RuntimeError("No valid events were loaded from the requested backgrounds.")

    all_values = np.concatenate(merged_values, axis=0)
    all_signed_weights = np.concatenate(merged_signed_weights, axis=0)
    optimization_weights = transform_optimization_weights(all_signed_weights, args.binning_weight_mode)

    edges, clipped_values, finite_mask = build_equal_content_edges(
        all_values,
        optimization_weights,
        args.n_bins,
        x_min=args.x_min,
        x_max=args.x_max,
    )

    clipped_signed_weights = all_signed_weights[finite_mask]
    clipped_optimization_weights = optimization_weights[finite_mask]
    bin_rows = histogram_summary(clipped_values, clipped_signed_weights, clipped_optimization_weights, edges)
    variable_config = make_variable_config(config_key, label, args.n_bins, edges)

    summary = {
        "sample_json": args.sample_json,
        "channels_requested": args.channels,
        "variable": args.variable,
        "variable_expression": variable_expression,
        "define_config": args.define_config,
        "config_key": config_key,
        "label": label,
        "n_bins": args.n_bins,
        "tree_name": args.tree_name,
        "weight_mode": args.weight_mode,
        "binning_weight_mode": args.binning_weight_mode,
        "lumi": lumi_value,
        "range": {
            "x_min": None if args.x_min is None else float(args.x_min),
            "x_max": None if args.x_max is None else float(args.x_max),
            "overflow_merged_to_last_bin": args.x_max is not None,
            "underflow_merged_to_first_bin": args.x_min is not None,
        },
        "bin_edges": [float(f"{edge:.12g}") for edge in edges],
        "channel_summaries": channel_weight_summaries,
        "bin_summaries": bin_rows,
        "variable_config": variable_config,
    }

    save_json(output_dir / args.summary_name, summary)
    save_json(output_dir / f"{config_key}_varConfig.json", variable_config)
    save_csv(output_dir / f"{config_key}_bin_summary.csv", bin_rows)

    print("optimized bin edges:", ", ".join(f"{edge:.12g}" for edge in edges))
    print(f"summary json: {output_dir / args.summary_name}")
    print(f"variable config json: {output_dir / f'{config_key}_varConfig.json'}")
    print(f"bin summary csv: {output_dir / f'{config_key}_bin_summary.csv'}")


if __name__ == "__main__":
    main()
