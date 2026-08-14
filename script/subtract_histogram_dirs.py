#!/usr/bin/env python3
import argparse
from pathlib import Path

import ROOT


ROOT.gROOT.SetBatch(True)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build output-channel histograms from (sum data channels over data dirs) - (sum MC channels over one MC dir)."
    )
    parser.add_argument("--data-dirs", nargs="+", required=True, help="Input data/DD histogram directories to add.")
    parser.add_argument("--mc-dir", required=True, help="Input MC histogram directory to subtract.")
    parser.add_argument("--output-dir", required=True, help="Output directory.")
    parser.add_argument("--data-channels", nargs="+", required=True, help="Channel names to sum from each data dir.")
    parser.add_argument("--mc-channels", nargs="+", required=True, help="Channel names to sum from the MC dir.")
    parser.add_argument("--output-channel", required=True, help="Output channel name prefix.")
    parser.add_argument(
        "--data-weight-prefix",
        default="DDTotalWeight",
        help="Nominal weight prefix used in data/DD files, e.g. DDTotalWeight.",
    )
    parser.add_argument(
        "--mc-weight-prefix",
        default="MCTotalWeight",
        help="Nominal weight prefix used in MC subtraction files, e.g. MCTotalWeight.",
    )
    parser.add_argument(
        "--variables",
        nargs="+",
        default=None,
        help="Optional list of variable names to keep, e.g. leadingMuon_pt subleadingMuon_pt.",
    )
    return parser.parse_args()


def open_hist(path, hist_name):
    root_file = ROOT.TFile.Open(str(path), "READ")
    if not root_file or root_file.IsZombie():
        return None
    hist = root_file.Get(hist_name)
    if not hist:
        root_file.Close()
        return None
    out = hist.Clone(hist_name + "_clone")
    out.SetDirectory(0)
    root_file.Close()
    return out


def add_hist(total, hist):
    if hist is None:
        return total
    if total is None:
        out = hist.Clone(hist.GetName() + "_sum")
        out.SetDirectory(0)
        return out
    total.Add(hist)
    return total


def collect_remainders(data_dirs, data_channels):
    remainders = set()
    for data_dir in data_dirs:
        directory = Path(data_dir)
        for channel in data_channels:
            prefix = channel + "_"
            for path in directory.glob(prefix + "*.root"):
                if not path.is_file():
                    continue
                stem = path.stem
                if not stem.startswith(prefix):
                    continue
                remainders.add(stem[len(prefix):])
    return sorted(remainders)


def filter_remainders_by_variables(remainders, variables):
    if not variables:
        return remainders
    variable_set = set(variables)
    selected = []
    for remainder in remainders:
        matched = False
        for variable in variable_set:
            if remainder == variable:
                matched = True
                break
            if remainder.startswith(variable + "_"):
                matched = True
                break
            if "_vs_" in remainder and remainder.startswith(variable + "_vs_"):
                matched = True
                break
        if matched:
            selected.append(remainder)
    return selected


def map_remainder_to_mc(remainder, data_weight_prefix, mc_weight_prefix):
    data_nominal = "_" + data_weight_prefix
    mc_nominal = "_" + mc_weight_prefix
    if remainder.endswith(data_nominal):
        return remainder[: -len(data_nominal)] + mc_nominal

    data_marker = data_nominal + "_"
    if data_marker in remainder:
        head, tail = remainder.rsplit(data_marker, 1)
        return head + mc_nominal + "_" + tail

    return remainder


def sum_channel_group(base_dirs, channels, remainder):
    combined = None
    used_files = []
    missing_files = []
    for base_dir in base_dirs:
        directory = Path(base_dir)
        for channel in channels:
            stem = channel + "_" + remainder
            path = directory / f"{stem}.root"
            hist = open_hist(path, stem)
            if hist is None:
                missing_files.append(str(path))
                continue
            used_files.append(str(path))
            combined = add_hist(combined, hist)
    return combined, used_files, missing_files


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    remainders = collect_remainders(args.data_dirs, args.data_channels)
    remainders = filter_remainders_by_variables(remainders, args.variables)
    if not remainders:
        raise RuntimeError("No input histograms were found for the requested data channels.")

    produced = 0
    for remainder in remainders:
        data_hist, used_data_files, missing_data_files = sum_channel_group(args.data_dirs, args.data_channels, remainder)
        if data_hist is None:
            continue

        mc_remainder = map_remainder_to_mc(remainder, args.data_weight_prefix, args.mc_weight_prefix)
        mc_hist, used_mc_files, missing_mc_files = sum_channel_group([args.mc_dir], args.mc_channels, mc_remainder)

        output_stem = args.output_channel + "_" + remainder
        output_hist = data_hist.Clone(output_stem)
        output_hist.SetDirectory(0)
        output_hist.SetTitle(output_stem)
        if mc_hist is not None:
            output_hist.Add(mc_hist, -1.0)

        output_path = output_dir / f"{output_stem}.root"
        output_file = ROOT.TFile.Open(str(output_path), "RECREATE")
        if not output_file or output_file.IsZombie():
            raise RuntimeError(f"Cannot create output ROOT file: {output_path}")
        output_hist.Write(output_stem)
        output_file.Close()

        produced += 1
        print(f"saved: {output_path}")
        print("  data inputs:")
        for file_path in used_data_files:
            print(f"    {file_path}")
        print("  mc inputs:")
        for file_path in used_mc_files:
            print(f"    {file_path}")
        if missing_data_files:
            print(f"  missing data files: {len(missing_data_files)}")
        if missing_mc_files:
            print(f"  missing mc files: {len(missing_mc_files)}")

    if produced == 0:
        raise RuntimeError("No output histograms were produced.")


if __name__ == "__main__":
    main()
