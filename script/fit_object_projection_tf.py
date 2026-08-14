#!/usr/bin/env python3
import argparse
import json
import math
from pathlib import Path

import ROOT


ROOT.gROOT.SetBatch(True)


def parse_bins(text):
    values = [float(item) for item in text.split(",") if item.strip()]
    if len(values) < 2:
        raise ValueError("binning must contain at least two edges")
    for left, right in zip(values, values[1:]):
        if right <= left:
            raise ValueError("bin edges must be strictly increasing")
    return values


def first_hist_name(root_file):
    for key in root_file.GetListOfKeys():
        obj = key.ReadObj()
        if obj.InheritsFrom("TH1"):
            return key.GetName()
    return None


def load_process_scale_map(path):
    if not path:
        return {}
    with Path(path).open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise RuntimeError("process scale json must be a dictionary")
    scale_map = {}
    for key, value in payload.items():
        try:
            scale_map[str(key)] = float(value)
        except (TypeError, ValueError) as exc:
            raise RuntimeError(f"Invalid process scale for {key}: {value}") from exc
    return scale_map


def load_hist(path):
    root_file = ROOT.TFile.Open(str(path), "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Cannot open ROOT file: {path}")
    hist_name = first_hist_name(root_file)
    if hist_name is None:
        root_file.Close()
        raise RuntimeError(f"No TH1 found in ROOT file: {path}")
    hist = root_file.Get(hist_name)
    out = hist.Clone(hist_name + "_clone")
    out.SetDirectory(0)
    root_file.Close()
    return out


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
    if suffix.startswith("vs_"):
        return False
    if weight is not None:
        return suffix == weight
    return suffix in {"one", "MCTotalWeight", "DDTotalWeight"}


def extract_process_name(stem, variable, weight, process=None):
    if process:
        return process
    if weight is not None:
        suffix = f"_{variable}_{weight}"
        if stem.endswith(suffix):
            return stem[: -len(suffix)]
    return stem.split("_", 1)[0]


def find_hist_files(directory, variable, process=None, weight=None):
    directory = Path(directory)
    if not directory.exists():
        raise RuntimeError(f"Input directory does not exist: {directory}")
    pattern = "*.root" if process is None else f"{process}_*.root"
    matches = [
        path for path in sorted(directory.glob(pattern))
        if variable_matches(path.stem, variable, process, weight)
    ]
    if not matches:
        raise RuntimeError(f"No histogram matched {directory} for variable={variable}, process={process}, weight={weight}")
    return matches


def same_1d_binning(a, b):
    if a.GetDimension() != 1 or b.GetDimension() != 1:
        return False
    if a.GetNbinsX() != b.GetNbinsX():
        return False
    axa = a.GetXaxis()
    axb = b.GetXaxis()
    for idx in range(1, axa.GetNbins() + 1):
        if abs(axa.GetBinLowEdge(idx) - axb.GetBinLowEdge(idx)) > 1e-9:
            return False
        if abs(axa.GetBinUpEdge(idx) - axb.GetBinUpEdge(idx)) > 1e-9:
            return False
    return True


def sum_from_dir(directory, variable, process=None, weight=None, return_components=False, process_scale_map=None):
    total = None
    components = []
    process_scale_map = process_scale_map or {}
    for path in find_hist_files(directory, variable, process, weight):
        hist = load_hist(path)
        if hist.GetDimension() != 1:
            raise RuntimeError(f"{path} is not a 1D histogram")
        process_name = extract_process_name(path.stem, variable, weight, process)
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
            if not same_1d_binning(total, hist):
                raise RuntimeError(f"Incompatible binning while summing {path}")
            total.Add(hist)
    if return_components:
        return total, components
    return total


def sum_processes_from_dir(directory, variable, processes, weight=None, return_components=False, process_scale_map=None):
    if not processes:
        raise RuntimeError("sum_processes_from_dir requires at least one process")
    total = None
    all_components = []
    for process in processes:
        hist, components = sum_from_dir(
            directory,
            variable,
            process=process,
            weight=weight,
            return_components=True,
            process_scale_map=process_scale_map,
        )
        all_components.extend(components)
        if total is None:
            total = hist.Clone(variable + "_selected_process_sum")
            total.SetDirectory(0)
        else:
            if not same_1d_binning(total, hist):
                raise RuntimeError(f"Incompatible binning while summing process {process} in {directory}")
            total.Add(hist)
    if return_components:
        return total, all_components
    return total


def sum_dirs(directories, variable, process=None, weight=None, return_components=False, process_scale_map=None):
    total = None
    all_components = []
    for directory in directories:
        if return_components:
            hist, components = sum_from_dir(directory, variable, process, weight, True, process_scale_map)
            all_components.extend(components)
        else:
            hist = sum_from_dir(directory, variable, process, weight, False, process_scale_map)
        if total is None:
            total = hist.Clone(variable + "_data_sum")
            total.SetDirectory(0)
        else:
            if not same_1d_binning(total, hist):
                raise RuntimeError(f"Incompatible binning while summing directory {directory}")
            total.Add(hist)
    if return_components:
        return total, all_components
    return total


def hist_bins(hist):
    bins = []
    axis = hist.GetXaxis()
    for ibin in range(1, hist.GetNbinsX() + 1):
        bins.append({
            "bin": ibin,
            "low": axis.GetBinLowEdge(ibin),
            "high": axis.GetBinUpEdge(ibin),
            "content": hist.GetBinContent(ibin),
            "error": hist.GetBinError(ibin),
        })
    return bins


def hist_summary(hist):
    return {
        "integral": hist.Integral(),
        "bins": hist_bins(hist),
    }


def component_summaries(components):
    out = []
    for path, hist in components:
        out.append({
            "path": path,
            "integral": hist.Integral(),
            "bins": hist_bins(hist),
        })
    return out


def build_residual(label, direct_dir, data_dirs, subtract_dir, variable, process, weight, subtract_process, subtract_processes, subtract_weight, diagnostics, process_scale_map=None):
    if direct_dir:
        direct, components = sum_from_dir(direct_dir, variable, process, weight, True, process_scale_map)
        diagnostics[label] = {
            "variable": variable,
            "mode": "direct",
            "direct": hist_summary(direct),
            "direct_components": component_summaries(components),
        }
        return direct
    if not data_dirs or not subtract_dir:
        raise RuntimeError(f"{label}: provide either direct-dir or data-dirs plus subtract-dir")

    data, data_components = sum_dirs(data_dirs, variable, process, weight, True, process_scale_map)
    if subtract_processes:
        subtract, subtract_components = sum_processes_from_dir(
            subtract_dir,
            variable,
            subtract_processes,
            subtract_weight or weight,
            True,
            process_scale_map,
        )
    else:
        subtract, subtract_components = sum_from_dir(
            subtract_dir,
            variable,
            subtract_process or process,
            subtract_weight or weight,
            True,
            process_scale_map,
        )
    if not same_1d_binning(data, subtract):
        raise RuntimeError(f"{label}: data and subtraction binning are incompatible")
    residual = data.Clone(label.replace(" ", "_") + "_residual")
    residual.SetDirectory(0)
    residual.Add(subtract, -1.0)
    diagnostics[label] = {
        "variable": variable,
        "mode": "data_minus_subtract",
        "data": hist_summary(data),
        "subtract": hist_summary(subtract),
        "residual": hist_summary(residual),
        "data_components": component_summaries(data_components),
        "subtract_components": component_summaries(subtract_components),
    }
    return residual


def find_tf_bin(tf_bins, x):
    for idx, (low, high) in enumerate(zip(tf_bins, tf_bins[1:])):
        if x >= low and x < high:
            return idx
    if abs(x - tf_bins[-1]) < 1e-9:
        return len(tf_bins) - 2
    return None


def invert_matrix(matrix):
    n = len(matrix)
    aug = []
    for i in range(n):
        aug.append([float(matrix[i][j]) for j in range(n)] + [1.0 if i == j else 0.0 for j in range(n)])

    for col in range(n):
        pivot = max(range(col, n), key=lambda row: abs(aug[row][col]))
        if abs(aug[pivot][col]) < 1e-14:
            raise RuntimeError("Singular normal matrix in TF fit")
        aug[col], aug[pivot] = aug[pivot], aug[col]
        scale = aug[col][col]
        for j in range(2 * n):
            aug[col][j] /= scale
        for row in range(n):
            if row == col:
                continue
            factor = aug[row][col]
            if factor == 0.0:
                continue
            for j in range(2 * n):
                aug[row][j] -= factor * aug[col][j]
    return [[aug[i][n + j] for j in range(n)] for i in range(n)]


def solve_weighted_ls(rows, values, sigmas):
    npar = len(rows[0])
    normal = [[0.0 for _ in range(npar)] for _ in range(npar)]
    rhs = [0.0 for _ in range(npar)]

    for row, value, sigma in zip(rows, values, sigmas):
        weight = 1.0 / (sigma * sigma) if sigma > 0 else 1.0
        for i in range(npar):
            rhs[i] += row[i] * value * weight
            for j in range(npar):
                normal[i][j] += row[i] * row[j] * weight

    cov = invert_matrix(normal)
    params = [sum(cov[i][j] * rhs[j] for j in range(npar)) for i in range(npar)]
    errors = [math.sqrt(max(cov[i][i], 0.0)) for i in range(npar)]
    return params, errors, cov


def add_hists(name, a, b):
    out = a.Clone(name)
    out.SetDirectory(0)
    out.Add(b)
    return out


def get_tf_average_weights(anti_hist, tf_bins):
    weights = [0.0 for _ in range(len(tf_bins) - 1)]
    total = 0.0
    for ibin in range(1, anti_hist.GetNbinsX() + 1):
        center = anti_hist.GetXaxis().GetBinCenter(ibin)
        tf_idx = find_tf_bin(tf_bins, center)
        if tf_idx is None:
            continue
        content = anti_hist.GetBinContent(ibin)
        weights[tf_idx] += content
        total += content
    if total <= 0:
        raise RuntimeError("Cannot compute TFmean: AR antiIso histogram has no yield inside TF bins")
    return [value / total for value in weights]


def weighted_tf_mean(params, average_weights):
    return sum(value * weight for value, weight in zip(params, average_weights))


def weighted_tf_mean_error(cov, average_weights):
    err2 = 0.0
    for i, wi in enumerate(average_weights):
        for j, wj in enumerate(average_weights):
            err2 += wi * cov[i][j] * wj
    return math.sqrt(max(err2, 0.0))


def make_prediction(name, template, anti_hist, tight_hist, tf_bins, params, average_weights):
    out = template.Clone(name)
    out.Reset("ICES")
    tf_mean = weighted_tf_mean(params, average_weights)
    for ibin in range(1, template.GetNbinsX() + 1):
        center = template.GetXaxis().GetBinCenter(ibin)
        tf_idx = find_tf_bin(tf_bins, center)
        anti_part = 0.0 if tf_idx is None else params[tf_idx] * anti_hist.GetBinContent(ibin)
        tight_part = tf_mean * tight_hist.GetBinContent(ibin)
        out.SetBinContent(ibin, anti_part + tight_part)
    return out


def build_system(sr_all, ar_anti, ar_tight, tf_bins, average_weights):
    if not same_1d_binning(sr_all, ar_anti) or not same_1d_binning(sr_all, ar_tight):
        raise RuntimeError("Input histograms must have identical 1D binning")

    n_tf = len(tf_bins) - 1
    rows = []
    values = []
    sigmas = []
    labels = []
    for ibin in range(1, sr_all.GetNbinsX() + 1):
        center = sr_all.GetXaxis().GetBinCenter(ibin)
        tf_idx = find_tf_bin(tf_bins, center)
        row = [ar_tight.GetBinContent(ibin) * weight for weight in average_weights]
        if tf_idx is not None:
            row[tf_idx] += ar_anti.GetBinContent(ibin)
        rows.append(row)
        values.append(sr_all.GetBinContent(ibin))
        sigma = sr_all.GetBinError(ibin)
        if sigma <= 0:
            sigma = math.sqrt(abs(sr_all.GetBinContent(ibin))) if sr_all.GetBinContent(ibin) > 0 else 1.0
        sigmas.append(sigma)
        labels.append({
            "bin": ibin,
            "center": center,
            "tf_bin": None if tf_idx is None else tf_idx + 1,
            "sr": sr_all.GetBinContent(ibin),
            "sigma": sigma,
            "anti": ar_anti.GetBinContent(ibin),
            "tight": ar_tight.GetBinContent(ibin),
            "tfmean_tight_coefficients": list(average_weights),
        })
    return rows, values, sigmas, labels


def build_region_system(sr_all, ar_anti, ar_tight, tf_bins, average_weights):
    if not same_1d_binning(sr_all, ar_anti) or not same_1d_binning(sr_all, ar_tight):
        raise RuntimeError("Input histograms must have identical 1D binning")

    n_tf = len(tf_bins) - 1
    anti_sums = [0.0 for _ in range(n_tf)]
    anti_err2 = [0.0 for _ in range(n_tf)]
    tight_sums = [0.0 for _ in range(n_tf)]
    tight_err2 = [0.0 for _ in range(n_tf)]
    sr_sums = [0.0 for _ in range(n_tf)]
    sr_err2 = [0.0 for _ in range(n_tf)]

    for ibin in range(1, sr_all.GetNbinsX() + 1):
        center = sr_all.GetXaxis().GetBinCenter(ibin)
        tf_idx = find_tf_bin(tf_bins, center)
        if tf_idx is None:
            continue
        anti_sums[tf_idx] += ar_anti.GetBinContent(ibin)
        anti_err2[tf_idx] += ar_anti.GetBinError(ibin) ** 2
        tight_sums[tf_idx] += ar_tight.GetBinContent(ibin)
        tight_err2[tf_idx] += ar_tight.GetBinError(ibin) ** 2
        sr_sums[tf_idx] += sr_all.GetBinContent(ibin)
        sr_err2[tf_idx] += sr_all.GetBinError(ibin) ** 2

    rows = []
    values = []
    sigmas = []
    labels = []
    sr_total = sum(sr_sums)
    anti_total = sum(anti_sums)
    tight_total = sum(tight_sums)
    for idx, (low, high) in enumerate(zip(tf_bins, tf_bins[1:])):
        row = [tight_sums[idx] * weight for weight in average_weights]
        row[idx] += anti_sums[idx]
        sigma = math.sqrt(sr_err2[idx]) if sr_err2[idx] > 0 else (math.sqrt(abs(sr_sums[idx])) if sr_sums[idx] > 0 else 1.0)
        rows.append(row)
        values.append(sr_sums[idx])
        sigmas.append(sigma)
        labels.append({
            "region": idx + 1,
            "low": low,
            "high": high,
            "tf_bin": idx + 1,
            "sr": sr_sums[idx],
            "sigma": sigma,
            "anti": anti_sums[idx],
            "tight": tight_sums[idx],
            "sr_fraction": 0.0 if sr_total <= 0 else sr_sums[idx] / sr_total,
            "anti_fraction": 0.0 if anti_total <= 0 else anti_sums[idx] / anti_total,
            "tight_fraction": 0.0 if tight_total <= 0 else tight_sums[idx] / tight_total,
            "tfmean_tight_coefficients": list(average_weights),
        })
    return rows, values, sigmas, labels


def parse_args():
    parser = argparse.ArgumentParser(
        description="Fit 1D object-projection TF: SR leading+subleading = TF(anti) * AR anti + <TF(anti)> * AR tight."
    )
    parser.add_argument("--ar-direct-dir", default=None)
    parser.add_argument("--ar-data-dirs", nargs="+", default=None)
    parser.add_argument("--ar-subtract-dir", default=None)
    parser.add_argument("--sr-direct-dir", default=None)
    parser.add_argument("--sr-data-dirs", nargs="+", default=None)
    parser.add_argument("--sr-subtract-dir", default=None)
    parser.add_argument("--anti-var", default="antiIsoMuon_pt")
    parser.add_argument("--tight-var", default="tightMuon_pt")
    parser.add_argument("--sr-leading-var", default="leadingMuon_pt")
    parser.add_argument("--sr-subleading-var", default="subleadingMuon_pt")
    parser.add_argument("--tf-bins", default="15,30,50,100,205")
    parser.add_argument("--process", default=None)
    parser.add_argument("--weight-suffix", default=None)
    parser.add_argument("--subtract-process", default=None)
    parser.add_argument("--subtract-processes", nargs="+", default=None)
    parser.add_argument("--subtract-weight-suffix", default=None)
    parser.add_argument("--process-scale-json", default=None, help="Optional json dictionary of per-process multiplicative scales.")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--tag", default="object_projection_tf")
    parser.add_argument("--fit-level", choices=["hist-bin", "tf-region"], default="hist-bin")
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()

    if bool(args.ar_data_dirs) != bool(args.ar_subtract_dir):
        parser.error("--ar-data-dirs and --ar-subtract-dir must be provided together")
    if bool(args.sr_data_dirs) != bool(args.sr_subtract_dir):
        parser.error("--sr-data-dirs and --sr-subtract-dir must be provided together")
    if bool(args.ar_direct_dir) == bool(args.ar_data_dirs):
        parser.error("provide exactly one AR input mode: --ar-direct-dir or --ar-data-dirs/--ar-subtract-dir")
    if bool(args.sr_direct_dir) == bool(args.sr_data_dirs):
        parser.error("provide exactly one SR input mode: --sr-direct-dir or --sr-data-dirs/--sr-subtract-dir")
    return args


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    tf_bins = parse_bins(args.tf_bins)
    diagnostics = {}
    process_scale_map = load_process_scale_map(args.process_scale_json)

    ar_anti = build_residual(
        "AR anti",
        args.ar_direct_dir,
        args.ar_data_dirs,
        args.ar_subtract_dir,
        args.anti_var,
        args.process,
        args.weight_suffix,
        args.subtract_process,
        args.subtract_processes,
        args.subtract_weight_suffix,
        diagnostics,
        process_scale_map,
    )
    ar_tight = build_residual(
        "AR tight",
        args.ar_direct_dir,
        args.ar_data_dirs,
        args.ar_subtract_dir,
        args.tight_var,
        args.process,
        args.weight_suffix,
        args.subtract_process,
        args.subtract_processes,
        args.subtract_weight_suffix,
        diagnostics,
        process_scale_map,
    )
    sr_leading = build_residual(
        "SR leading",
        args.sr_direct_dir,
        args.sr_data_dirs,
        args.sr_subtract_dir,
        args.sr_leading_var,
        args.process,
        args.weight_suffix,
        args.subtract_process,
        args.subtract_processes,
        args.subtract_weight_suffix,
        diagnostics,
        process_scale_map,
    )
    sr_subleading = build_residual(
        "SR subleading",
        args.sr_direct_dir,
        args.sr_data_dirs,
        args.sr_subtract_dir,
        args.sr_subleading_var,
        args.process,
        args.weight_suffix,
        args.subtract_process,
        args.subtract_processes,
        args.subtract_weight_suffix,
        diagnostics,
        process_scale_map,
    )
    sr_all = add_hists(args.tag + "_sr_all", sr_leading, sr_subleading)

    average_weights = get_tf_average_weights(ar_anti, tf_bins)
    if args.fit_level == "hist-bin":
        rows, values, sigmas, labels = build_system(sr_all, ar_anti, ar_tight, tf_bins, average_weights)
    else:
        rows, values, sigmas, labels = build_region_system(sr_all, ar_anti, ar_tight, tf_bins, average_weights)
    params, errors, cov = solve_weighted_ls(rows, values, sigmas)
    prediction = make_prediction(args.tag + "_prediction", sr_all, ar_anti, ar_tight, tf_bins, params, average_weights)

    chi2 = 0.0
    if args.fit_level == "hist-bin":
        for ibin in range(1, sr_all.GetNbinsX() + 1):
            sigma = sigmas[ibin - 1]
            diff = prediction.GetBinContent(ibin) - sr_all.GetBinContent(ibin)
            chi2 += diff * diff / (sigma * sigma) if sigma > 0 else 0.0
    else:
        for entry, row, sigma in zip(labels, rows, sigmas):
            pred = sum(coeff * param for coeff, param in zip(row, params))
            diff = pred - entry["sr"]
            chi2 += diff * diff / (sigma * sigma) if sigma > 0 else 0.0
    ndof = len(values) - len(params)

    anti_weighted_tf_mean = weighted_tf_mean(params, average_weights)
    anti_weighted_tf_mean_error = weighted_tf_mean_error(cov, average_weights)

    result = {
        "tf_bins": tf_bins,
        "tf_values": params,
        "tf_errors": errors,
        "tfmean_weights": average_weights,
        "anti_weighted_tf_mean": anti_weighted_tf_mean,
        "anti_weighted_tf_mean_error": anti_weighted_tf_mean_error,
        "fit_level": args.fit_level,
        "process_scale_json": args.process_scale_json,
        "process_scales": process_scale_map,
        "chi2": chi2,
        "ndof": ndof,
        "equations": [],
        "covariance": cov,
        "inputs": diagnostics,
    }
    for label, row in zip(labels, rows):
        pred = sum(coeff * param for coeff, param in zip(row, params))
        entry = dict(label)
        entry["prediction"] = pred
        entry["pull"] = (pred - entry["sr"]) / entry["sigma"] if entry["sigma"] > 0 else 0.0
        result["equations"].append(entry)

    json_path = output_dir / f"{args.tag}.json"
    with open(json_path, "w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)

    root_path = output_dir / f"{args.tag}.root"
    fout = ROOT.TFile.Open(str(root_path), "RECREATE")
    ar_anti.Clone(args.tag + "_ar_anti").Write(args.tag + "_ar_anti")
    ar_tight.Clone(args.tag + "_ar_tight").Write(args.tag + "_ar_tight")
    sr_leading.Clone(args.tag + "_sr_leading").Write(args.tag + "_sr_leading")
    sr_subleading.Clone(args.tag + "_sr_subleading").Write(args.tag + "_sr_subleading")
    sr_all.Clone(args.tag + "_sr_all").Write(args.tag + "_sr_all")
    prediction.Write(args.tag + "_prediction")
    tf_hist = ROOT.TH1D(args.tag + "_tf", args.tag + "_tf", len(tf_bins) - 1, array_to_c(tf_bins))
    for idx, (value, error) in enumerate(zip(params, errors), start=1):
        tf_hist.SetBinContent(idx, value)
        tf_hist.SetBinError(idx, error)
    tf_hist.Write()
    fout.Close()

    print("Fit result")
    print(f"  fit level: {args.fit_level}")
    print(f"  chi2/ndof: {chi2:.6f} / {ndof}")
    for idx, (low, high, value, error) in enumerate(zip(tf_bins, tf_bins[1:], params, errors), start=1):
        print(f"  TF bin {idx} [{low}, {high}): {value:.8f} +- {error:.8f}")
    print(f"  anti-weighted TFmean: {anti_weighted_tf_mean:.8f} +- {anti_weighted_tf_mean_error:.8f}")
    print(f"Wrote {json_path}")
    print(f"Wrote {root_path}")

    if args.debug:
        print("Input diagnostics")
        for label, info in diagnostics.items():
            print(f"  {label} ({info['variable']}):")
            if info["mode"] == "direct":
                print(f"    direct integral={info['direct']['integral']:.6f}")
            else:
                print(f"    data integral={info['data']['integral']:.6f}")
                print(f"    subtract integral={info['subtract']['integral']:.6f}")
                print(f"    residual integral={info['residual']['integral']:.6f}")
                for bin_info in info["residual"]["bins"]:
                    if bin_info["content"] < 0:
                        print(
                            f"    NEG residual bin {bin_info['bin']} "
                            f"[{bin_info['low']:.3f}, {bin_info['high']:.3f}): "
                            f"{bin_info['content']:.6f} +- {bin_info['error']:.6f}"
                        )
            print("    subtract components:")
            for component in info.get("subtract_components", []):
                print(f"      {Path(component['path']).name}: integral={component['integral']:.6f}")

        print("Equation debug")
        for entry in result["equations"]:
            if args.fit_level == "hist-bin":
                print(
                    f"  bin {entry['bin']} center={entry['center']:.3f}: "
                    f"obs={entry['sr']:.6f}, pred={entry['prediction']:.6f}, "
                    f"pull={entry['pull']:.3f}, tf_bin={entry['tf_bin']}"
                )
            else:
                print(
                    f"  region {entry['region']} [{entry['low']:.3f}, {entry['high']:.3f}): "
                    f"obs={entry['sr']:.6f}, pred={entry['prediction']:.6f}, pull={entry['pull']:.3f}, "
                    f"anti={entry['anti']:.6f} ({100.0 * entry['anti_fraction']:.2f}%), "
                    f"tight={entry['tight']:.6f} ({100.0 * entry['tight_fraction']:.2f}%), "
                    f"sr={entry['sr']:.6f} ({100.0 * entry['sr_fraction']:.2f}%)"
                )


def array_to_c(values):
    from array import array
    return array("d", values)


if __name__ == "__main__":
    main()
