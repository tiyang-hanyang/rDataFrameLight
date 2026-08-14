#!/usr/bin/env python3
"""Fit a two-dimensional object-projection transfer factor.

The fitted model is

  SR(x, y) = TF(x, y) * AR_anti(x, y) + <TF> * AR_tight(x, y),

where <TF> is the AR antiIso-yield-weighted average of the coarse 2D TF.
"""

import argparse
import json
import math
import sys
from array import array
from pathlib import Path

import ROOT


ROOT.gROOT.SetBatch(True)

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import fit_object_projection_tf as one_d  # noqa: E402


def parse_bins(text):
    values = [float(item.strip()) for item in text.split(",") if item.strip()]
    if len(values) < 2 or any(right <= left for left, right in zip(values, values[1:])):
        raise ValueError("binning must contain strictly increasing edges")
    return values


def variable_matches(stem, variable, process, weight):
    if process:
        prefix = process + "_"
        if not stem.startswith(prefix):
            return False
        remainder = stem[len(prefix):]
    else:
        split_pos = stem.find("_")
        if split_pos < 0:
            return False
        remainder = stem[split_pos + 1:]
    prefix = variable + "_"
    if not remainder.startswith(prefix):
        return False
    suffix = remainder[len(prefix):]
    if weight is not None:
        return suffix == weight
    return suffix in {"one", "MCTotalWeight", "DDTotalWeight"}


def same_2d_binning(a, b):
    if a.GetDimension() != 2 or b.GetDimension() != 2:
        return False
    if a.GetNbinsX() != b.GetNbinsX() or a.GetNbinsY() != b.GetNbinsY():
        return False
    for axis_a, axis_b in ((a.GetXaxis(), b.GetXaxis()), (a.GetYaxis(), b.GetYaxis())):
        for index in range(1, axis_a.GetNbins() + 1):
            if abs(axis_a.GetBinLowEdge(index) - axis_b.GetBinLowEdge(index)) > 1e-9:
                return False
            if abs(axis_a.GetBinUpEdge(index) - axis_b.GetBinUpEdge(index)) > 1e-9:
                return False
    return True


# Reuse the established input/subtraction implementation, but allow 2D names
# and require compatible 2D histograms instead of 1D histograms.
one_d.variable_matches = variable_matches
one_d.same_1d_binning = same_2d_binning


def sum_from_dir_2d(directory, variable, process=None, weight=None, return_components=False, process_scale_map=None):
    total = None
    components = []
    process_scale_map = process_scale_map or {}
    for path in one_d.find_hist_files(directory, variable, process, weight):
        hist = one_d.load_hist(path)
        if hist.GetDimension() != 2:
            raise RuntimeError(f"{path} is not a 2D histogram")
        process_name = one_d.extract_process_name(path.stem, variable, weight, process)
        scale = process_scale_map.get(process_name, 1.0)
        if scale != 1.0:
            hist.Scale(scale)
        component = hist.Clone(path.stem)
        component.SetDirectory(0)
        components.append((str(path), component))
        if total is None:
            total = hist.Clone(hist.GetName() + "_sum")
            total.SetDirectory(0)
        else:
            if not same_2d_binning(total, hist):
                raise RuntimeError(f"Incompatible 2D binning while summing {path}")
            total.Add(hist)
    if return_components:
        return total, components
    return total


def sum_processes_from_dir_2d(directory, variable, processes, weight=None, return_components=False, process_scale_map=None):
    total = None
    components = []
    for process in processes:
        hist, process_components = sum_from_dir_2d(
            directory, variable, process, weight, True, process_scale_map
        )
        components.extend(process_components)
        if total is None:
            total = hist.Clone(variable + "_selected_process_sum")
            total.SetDirectory(0)
        else:
            if not same_2d_binning(total, hist):
                raise RuntimeError(f"Incompatible 2D binning while summing process {process}")
            total.Add(hist)
    if return_components:
        return total, components
    return total


def sum_dirs_2d(directories, variable, process=None, weight=None, return_components=False, process_scale_map=None):
    total = None
    components = []
    for directory in directories:
        hist, directory_components = sum_from_dir_2d(
            directory, variable, process, weight, True, process_scale_map
        )
        components.extend(directory_components)
        if total is None:
            total = hist.Clone(variable + "_data_sum")
            total.SetDirectory(0)
        else:
            if not same_2d_binning(total, hist):
                raise RuntimeError(f"Incompatible 2D binning while summing {directory}")
            total.Add(hist)
    if return_components:
        return total, components
    return total


one_d.sum_from_dir = sum_from_dir_2d
one_d.sum_processes_from_dir = sum_processes_from_dir_2d
one_d.sum_dirs = sum_dirs_2d


def find_tf_bin(edges, value):
    for index, (low, high) in enumerate(zip(edges, edges[1:])):
        if low <= value < high:
            return index
    if abs(value - edges[-1]) < 1e-9:
        return len(edges) - 2
    return None


def cell_index(x_edges, y_edges, x, y):
    ix = find_tf_bin(x_edges, x)
    iy = find_tf_bin(y_edges, y)
    if ix is None or iy is None:
        return None
    return iy * (len(x_edges) - 1) + ix


def average_weights(anti, x_edges, y_edges):
    n_x = len(x_edges) - 1
    weights = [0.0] * (n_x * (len(y_edges) - 1))
    total = 0.0
    for ix in range(1, anti.GetNbinsX() + 1):
        x = anti.GetXaxis().GetBinCenter(ix)
        for iy in range(1, anti.GetNbinsY() + 1):
            y = anti.GetYaxis().GetBinCenter(iy)
            index = cell_index(x_edges, y_edges, x, y)
            if index is None:
                continue
            content = anti.GetBinContent(ix, iy)
            weights[index] += content
            total += content
    if total <= 0:
        raise RuntimeError("AR antiIso histogram has no yield inside TF bins")
    return [value / total for value in weights]


def build_system(sr, anti, tight, x_edges, y_edges, weights):
    if not same_2d_binning(sr, anti) or not same_2d_binning(sr, tight):
        raise RuntimeError("Input histograms must have identical 2D binning")
    n_params = len(weights)
    rows, values, sigmas, labels = [], [], [], []
    for ix in range(1, sr.GetNbinsX() + 1):
        x = sr.GetXaxis().GetBinCenter(ix)
        for iy in range(1, sr.GetNbinsY() + 1):
            y = sr.GetYaxis().GetBinCenter(iy)
            index = cell_index(x_edges, y_edges, x, y)
            row = [tight.GetBinContent(ix, iy) * weight for weight in weights]
            if index is not None:
                row[index] += anti.GetBinContent(ix, iy)
            value = sr.GetBinContent(ix, iy)
            sigma = sr.GetBinError(ix, iy)
            if sigma <= 0:
                if value == 0:
                    # An empty residual bin with no uncertainty carries no
                    # usable constraint and must not be treated as sigma=1.
                    continue
                sigma = math.sqrt(abs(value))
            rows.append(row)
            values.append(value)
            sigmas.append(sigma)
            labels.append({
                "x_bin": ix,
                "y_bin": iy,
                "x_center": x,
                "y_center": y,
                "tf_bin": None if index is None else index + 1,
                "sr": value,
                "sigma": sigma,
                "anti": anti.GetBinContent(ix, iy),
                "tight": tight.GetBinContent(ix, iy),
            })
    return rows, values, sigmas, labels


def make_prediction(name, template, anti, tight, x_edges, y_edges, params, weights):
    out = template.Clone(name)
    out.Reset("ICES")
    tf_mean = one_d.weighted_tf_mean(params, weights)
    for ix in range(1, template.GetNbinsX() + 1):
        x = template.GetXaxis().GetBinCenter(ix)
        for iy in range(1, template.GetNbinsY() + 1):
            y = template.GetYaxis().GetBinCenter(iy)
            index = cell_index(x_edges, y_edges, x, y)
            anti_part = 0.0 if index is None else params[index] * anti.GetBinContent(ix, iy)
            tight_part = tf_mean * tight.GetBinContent(ix, iy)
            out.SetBinContent(ix, iy, anti_part + tight_part)
    return out, tf_mean


def summarize_tf_regions(hists, x_edges, y_edges):
    """Return yields and propagated bin-statistical errors per coarse TF cell."""
    summaries = []
    n_x = len(x_edges) - 1
    for iy, (y_low, y_high) in enumerate(zip(y_edges, y_edges[1:])):
        for ix, (x_low, x_high) in enumerate(zip(x_edges, x_edges[1:])):
            values = {}
            for name, hist in hists.items():
                total = 0.0
                error2 = 0.0
                for bin_x in range(1, hist.GetNbinsX() + 1):
                    x = hist.GetXaxis().GetBinCenter(bin_x)
                    if not (x_low <= x < x_high):
                        continue
                    for bin_y in range(1, hist.GetNbinsY() + 1):
                        y = hist.GetYaxis().GetBinCenter(bin_y)
                        if not (y_low <= y < y_high):
                            continue
                        total += hist.GetBinContent(bin_x, bin_y)
                        error2 += hist.GetBinError(bin_x, bin_y) ** 2
                values[name] = {"yield": total, "error": math.sqrt(max(error2, 0.0))}
            summaries.append({
                "tf_bin": iy * n_x + ix + 1,
                "x_range": [x_low, x_high],
                "y_range": [y_low, y_high],
                "values": values,
            })
    return summaries


def aggregate_tf_regions(rows, values, sigmas, labels, n_params):
    grouped = {}
    for row, value, sigma, label in zip(rows, values, sigmas, labels):
        index = label["tf_bin"]
        if index is None:
            continue
        slot = grouped.setdefault(index - 1, {
            "row": [0.0] * n_params,
            "value": 0.0,
            "sigma2": 0.0,
            "tf_bin": index,
        })
        slot["row"] = [left + right for left, right in zip(slot["row"], row)]
        slot["value"] += value
        slot["sigma2"] += sigma * sigma
    ordered = [grouped[index] for index in sorted(grouped)]
    out_rows = [entry["row"] for entry in ordered]
    out_values = [entry["value"] for entry in ordered]
    out_sigmas = [math.sqrt(entry["sigma2"]) for entry in ordered]
    out_labels = [{"tf_bin": entry["tf_bin"], "sr": entry["value"], "sigma": sigma}
                  for entry, sigma in zip(ordered, out_sigmas)]
    return out_rows, out_values, out_sigmas, out_labels


def parse_args():
    parser = argparse.ArgumentParser(description="Fit a 2D object-projection TF.")
    parser.add_argument("--ar-direct-dir", default=None)
    parser.add_argument("--ar-data-dirs", nargs="+", default=None)
    parser.add_argument("--ar-subtract-dir", default=None)
    parser.add_argument("--sr-direct-dir", default=None)
    parser.add_argument("--sr-data-dirs", nargs="+", default=None)
    parser.add_argument("--sr-subtract-dir", default=None)
    parser.add_argument("--anti-var", default="antiIsoMuon_pt_vs_antiIsoMuon_eta_2bin")
    parser.add_argument("--tight-var", default="tightMuon_pt_vs_tightMuon_eta_2bin")
    parser.add_argument("--sr-leading-var", default="leadingMuon_pt_vs_leadingMuon_eta_2bin")
    parser.add_argument("--sr-subleading-var", default="subleadingMuon_pt_vs_subleadingMuon_eta_2bin")
    parser.add_argument("--tf-bins-x", required=True)
    parser.add_argument("--tf-bins-y", required=True)
    parser.add_argument("--process", default=None)
    parser.add_argument("--weight-suffix", default=None)
    parser.add_argument("--subtract-process", default=None)
    parser.add_argument("--subtract-processes", nargs="+", default=None)
    parser.add_argument("--subtract-weight-suffix", default=None)
    parser.add_argument("--process-scale-json", default=None)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--tag", default="object_projection_tf_2d")
    parser.add_argument(
        "--fit-level", choices=["hist-bin", "tf-region"], default="tf-region",
        help="Fit coarse 2D yields by default; use hist-bin only as a fine-shape diagnostic.",
    )
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()
    if bool(args.ar_data_dirs) != bool(args.ar_subtract_dir):
        parser.error("--ar-data-dirs and --ar-subtract-dir must be provided together")
    if bool(args.sr_data_dirs) != bool(args.sr_subtract_dir):
        parser.error("--sr-data-dirs and --sr-subtract-dir must be provided together")
    if bool(args.ar_direct_dir) == bool(args.ar_data_dirs):
        parser.error("provide exactly one AR input mode")
    if bool(args.sr_direct_dir) == bool(args.sr_data_dirs):
        parser.error("provide exactly one SR input mode")
    return args


def main():
    args = parse_args()
    x_edges = parse_bins(args.tf_bins_x)
    y_edges = parse_bins(args.tf_bins_y)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    diagnostics = {}
    scales = one_d.load_process_scale_map(args.process_scale_json)

    def residual(label, direct, data, subtract):
        return one_d.build_residual(
            label, direct, data, subtract, args.anti_var if label == "AR anti" else args.tight_var,
            args.process, args.weight_suffix, args.subtract_process, args.subtract_processes,
            args.subtract_weight_suffix, diagnostics, scales,
        )

    ar_anti = residual("AR anti", args.ar_direct_dir, args.ar_data_dirs, args.ar_subtract_dir)
    ar_tight = one_d.build_residual(
        "AR tight", args.ar_direct_dir, args.ar_data_dirs, args.ar_subtract_dir,
        args.tight_var, args.process, args.weight_suffix, args.subtract_process,
        args.subtract_processes, args.subtract_weight_suffix, diagnostics, scales,
    )
    sr_leading = one_d.build_residual(
        "SR leading", args.sr_direct_dir, args.sr_data_dirs, args.sr_subtract_dir,
        args.sr_leading_var, args.process, args.weight_suffix, args.subtract_process,
        args.subtract_processes, args.subtract_weight_suffix, diagnostics, scales,
    )
    sr_subleading = one_d.build_residual(
        "SR subleading", args.sr_direct_dir, args.sr_data_dirs, args.sr_subtract_dir,
        args.sr_subleading_var, args.process, args.weight_suffix, args.subtract_process,
        args.subtract_processes, args.subtract_weight_suffix, diagnostics, scales,
    )
    sr_all = one_d.add_hists("sr_all", sr_leading, sr_subleading)
    weights = average_weights(ar_anti, x_edges, y_edges)
    fine_rows, fine_values, fine_sigmas, fine_labels = build_system(
        sr_all, ar_anti, ar_tight, x_edges, y_edges, weights
    )
    rows, values, sigmas, labels = fine_rows, fine_values, fine_sigmas, fine_labels
    if args.fit_level == "tf-region":
        rows, values, sigmas, labels = aggregate_tf_regions(
            rows, values, sigmas, labels, len(weights)
        )
    params, errors, covariance = one_d.solve_weighted_ls(rows, values, sigmas)
    prediction, tf_mean = make_prediction("prediction", sr_all, ar_anti, ar_tight, x_edges, y_edges, params, weights)
    chi2 = sum(((pred - value) / sigma) ** 2 for pred, value, sigma in zip(
        [sum(c * p for c, p in zip(row, params)) for row in rows], values, sigmas
    ) if sigma > 0)
    ndof = len(values) - len(params)
    fine_predictions = [sum(c * p for c, p in zip(row, params)) for row in fine_rows]
    fine_chi2 = sum(
        ((pred - value) / sigma) ** 2
        for pred, value, sigma in zip(fine_predictions, fine_values, fine_sigmas)
        if sigma > 0
    )
    fine_ndof = len(fine_values) - len(params)
    tfmean_error = one_d.weighted_tf_mean_error(covariance, weights)
    region_summaries = summarize_tf_regions({
        "ar_anti": ar_anti,
        "ar_tight": ar_tight,
        "sr_leading": sr_leading,
        "sr_subleading": sr_subleading,
        "sr_all": sr_all,
        "prediction": prediction,
    }, x_edges, y_edges)

    for label, row in zip(labels, rows):
        label["prediction"] = sum(c * p for c, p in zip(row, params))
        label["pull"] = (label["prediction"] - label["sr"]) / label["sigma"] if label["sigma"] > 0 else 0.0
        label["ratio"] = label["prediction"] / label["sr"] if label["sr"] != 0 else None

    result = {
        "tf_bins_x": x_edges,
        "tf_bins_y": y_edges,
        "tf_values": params,
        "tf_errors": errors,
        "tfmean_weights": weights,
        "anti_weighted_tf_mean": tf_mean,
        "anti_weighted_tf_mean_error": tfmean_error,
        "fit_level": args.fit_level,
        "chi2": chi2,
        "ndof": ndof,
        "fine_chi2": fine_chi2,
        "fine_ndof": fine_ndof,
        "equations": labels,
        "covariance": covariance,
        "region_summaries": region_summaries,
        "inputs": diagnostics,
    }
    json_path = output_dir / f"{args.tag}.json"
    json_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    root_path = output_dir / f"{args.tag}.root"
    fout = ROOT.TFile.Open(str(root_path), "RECREATE")
    for name, hist in (("ar_anti", ar_anti), ("ar_tight", ar_tight), ("sr_leading", sr_leading),
                       ("sr_subleading", sr_subleading), ("sr_all", sr_all), ("prediction", prediction)):
        hist.Write(name)
    tf_hist = ROOT.TH2D(
        "tf", "tf", len(x_edges) - 1, array("d", x_edges),
        len(y_edges) - 1, array("d", y_edges)
    )
    for iy in range(len(y_edges) - 1):
        for ix in range(len(x_edges) - 1):
            index = iy * (len(x_edges) - 1) + ix
            tf_hist.SetBinContent(ix + 1, iy + 1, params[index])
            tf_hist.SetBinError(ix + 1, iy + 1, errors[index])
    tf_hist.Write()
    fout.Close()

    print("Fit result")
    print(f"  fit level: {args.fit_level}")
    print(f"  chi2/ndof: {chi2:.6f} / {ndof}")
    if args.fit_level == "tf-region":
        print(f"  fine-bin closure chi2/ndof: {fine_chi2:.6f} / {fine_ndof}")
    for iy in range(len(y_edges) - 1):
        for ix in range(len(x_edges) - 1):
            index = iy * (len(x_edges) - 1) + ix
            print(f"  TF bin ({ix + 1}, {iy + 1}) [{x_edges[ix]}, {x_edges[ix + 1]}) x [{y_edges[iy]}, {y_edges[iy + 1]}): {params[index]:.8f} +- {errors[index]:.8f}")
    print(f"  anti-weighted TFmean: {tf_mean:.8f} +- {tfmean_error:.8f}")
    print("  coarse-region yields (yield +- stat. error):")
    for region in region_summaries:
        values = region["values"]
        print(
            f"    TF bin {region['tf_bin']} "
            f"x=[{region['x_range'][0]}, {region['x_range'][1]}) "
            f"y=[{region['y_range'][0]}, {region['y_range'][1]}): "
            + ", ".join(
                f"{name}={item['yield']:.6f} +- {item['error']:.6f}"
                for name, item in values.items()
            )
        )
    print(f"Wrote {json_path}")
    print(f"Wrote {root_path}")


if __name__ == "__main__":
    main()
