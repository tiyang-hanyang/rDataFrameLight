#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

import ROOT


ROOT.gROOT.SetBatch(True)

DEFAULT_LUMI_CONFIG = Path(__file__).resolve().parent.parent / "json" / "Lumi" / "Run3.json"
DEFAULT_XS_SYST_CONFIG = Path(__file__).resolve().parent.parent / "json" / "syst" / "Run3_XSonly.json"
DEFAULT_SYST_ALIAS_CONFIG = Path(__file__).resolve().parent.parent / "json" / "general_config" / "syst_NP.json"
DEFAULT_SHAPE_SYST_GROUP_CONFIG = (
    Path(__file__).resolve().parent.parent / "json" / "general_config" / "shape_syst_groups.json"
)


CATEGORIES = ["signal", "ttX", "ttb", "TTTT"]
PROCESSES = ["TTHH", "TTW", "TTH", "TTZ", "TTTT", "ttbar", "TTBB", "TV", "VV", "VVV"]
DEFAULT_BACKGROUND_PROCESSES = PROCESSES[1:]
B_TAG_SYSTS = [
    "btag_lf",
    "btag_lfstats1",
    "btag_lfstats2",
    "btag_hf",
    "btag_hfstats1",
    "btag_hfstats2",
    "btag_cferr1",
    "btag_cferr2",
    "btag_bfragmentation",
    "btag_fsrdef",
    "btag_hdamp",
    "btag_isrdef",
    "btag_jer",
    "btag_jes",
    "btag_muf",
    "btag_mur",
    "btag_pdfas",
    "btag_pileup",
    "btag_statistic",
    "btag_topmass",
    "btag_type3",
]

DEFAULT_FULL_SHAPE_SYST_GROUP = "all_corr_fullRun"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge signal_bestScore histograms across campaigns and generate a combined AsymptoticLimits card."
    )
    parser.add_argument(
        "--campaigns",
        nargs="+",
        required=True,
        help="Campaign labels, e.g. Run3Summer22NanoAODv12 Run3Summer22EENanoAODv12 Run3Summer23NanoAODv12 Run3Summer23BPixNanoAODv12 RunIII2024Summer24NanoAODv15",
    )
    parser.add_argument(
        "--input-dir",
        default=".",
        help="Base directory containing S4_dimuon_DNN_*_region_<campaign> folders.",
    )
    parser.add_argument(
        "--output-tag",
        default="Run2x3x4Combined",
        help="Tag appended to merged output folders and card names.",
    )
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Directory where merged folders and datacard are written.",
    )
    parser.add_argument(
        "--signal-campaigns",
        nargs="+",
        default=None,
        help=(
            "Optional campaign labels used only for the signal process TTHH. "
            "Use this when backgrounds follow --campaigns but signal exists only in a different campaign."
        ),
    )
    parser.add_argument(
        "--background-processes",
        nargs="+",
        default=list(DEFAULT_BACKGROUND_PROCESSES),
        metavar="PROCESS",
        help=(
            "Background channels to include. Only TTHH is treated as the signal; "
            "all other channels are read exclusively from this list. This may "
            "include custom DD channels such as nonprompt_2D_shapeTF_final."
        ),
    )
    parser.add_argument(
        "--dd-processes",
        nargs="+",
        default=[],
        metavar="PROCESS",
        help=(
            "Selected background channels whose input histograms use the "
            "DDTotalWeight suffix instead of MCTotalWeight. Every DD process "
            "must also be listed in --background-processes."
        ),
    )
    parser.add_argument(
        "--lumi-config",
        default=str(DEFAULT_LUMI_CONFIG),
        help="Path to lumi json used to derive campaign lumi rescaling.",
    )
    parser.add_argument(
        "--lumi-lnN",
        default="1.10",
        help="Flat lumi lnN value written to the card.",
    )
    parser.add_argument(
        "--xs-syst-config",
        default=str(DEFAULT_XS_SYST_CONFIG),
        help=(
            "Path to the process-wise XS lnN json. The file should be keyed by merged datacard "
            "process names such as TTBB, TTW, TTZ, ttbar."
        ),
    )
    parser.add_argument(
        "--disable-xs-lnN",
        action="store_true",
        help="Do not write process-wise XS lnN nuisances.",
    )
    parser.add_argument(
        "--categories",
        nargs="+",
        choices=CATEGORIES,
        default=list(CATEGORIES),
        help="Subset of DNN regions to include in the fit. Default: all four regions.",
    )
    parser.add_argument(
        "--flat-input",
        action="store_true",
        help="Read templates directly from a flat input directory instead of per-category subdirectories.",
    )
    parser.add_argument(
        "--single-region-bin",
        default=None,
        help="Datacard bin name for flat single-region input, e.g. SR or optimize.",
    )
    parser.add_argument(
        "--hist-base-name",
        default="Signal_bestScore",
        help="Histogram base name used in template ROOT names, e.g. Signal_bestScore or optimize_score.",
    )
    parser.add_argument(
        "--folder-template",
        default=None,
        help=(
            "Optional input folder template for non-flat inputs. Available fields: {category}, {campaign}. "
            "Example: S4_dimuon_DNN_mixed_{campaign}."
        ),
    )
    parser.add_argument(
        "--stat-only",
        action="store_true",
        help="Write a statistical-only card: no lumi lnN, no XS lnN, and no shape systematics.",
    )
    parser.add_argument(
        "--include-shape-syst",
        action="store_true",
        help="Include generic shape systematics from MCTotalWeight_<syst>_up/down histograms.",
    )
    parser.add_argument(
        "--shape-systs",
        nargs="+",
        default=None,
        help="Generic shape systematic names to include, e.g. CMS_scale_j_Absolute CMS_scale_j_Absolute_2024.",
    )
    parser.add_argument(
        "--shape-syst-groups",
        nargs="+",
        default=None,
        help=(
            "Named shape systematic groups defined in the shape group json. "
            f"Default for --include-shape-syst: {DEFAULT_FULL_SHAPE_SYST_GROUP}."
        ),
    )
    parser.add_argument(
        "--shape-syst-alias-config",
        default=str(DEFAULT_SYST_ALIAS_CONFIG),
        help="Path to the systematic alias json, e.g. json/general_config/syst_NP.json.",
    )
    parser.add_argument(
        "--shape-syst-group-config",
        default=str(DEFAULT_SHAPE_SYST_GROUP_CONFIG),
        help="Path to the shape systematic group json.",
    )
    parser.add_argument(
        "--include-btag-syst",
        action="store_true",
        help="Backward-compatible alias of --include-shape-syst for b-tag weight systematics.",
    )
    parser.add_argument(
        "--btag-systs",
        nargs="+",
        default=list(B_TAG_SYSTS),
        help="Backward-compatible alias list for --shape-systs.",
    )
    parser.add_argument(
        "--overwrite-output",
        action="store_true",
        help="Allow overwriting existing output ROOT files, merged directories, and datacard.",
    )
    return parser.parse_args()


def category_dir(base_dir, category, campaign, folder_template=None):
    if folder_template:
        return base_dir / folder_template.format(category=category, campaign=campaign)
    if category == "TTTT":
        name = f"S4_dimuon_DNN_TTTT_region_{campaign}"
    else:
        name = f"S4_dimuon_DNN_{category}_region_{campaign}"
    return base_dir / name


def merge_underflow_overflow(hist):
    nbins = hist.GetNbinsX()

    underflow_content = hist.GetBinContent(0)
    underflow_error = hist.GetBinError(0)
    first_content = hist.GetBinContent(1)
    first_error = hist.GetBinError(1)
    hist.SetBinContent(1, first_content + underflow_content)
    hist.SetBinError(1, (first_error * first_error + underflow_error * underflow_error) ** 0.5)
    hist.SetBinContent(0, 0.0)
    hist.SetBinError(0, 0.0)

    overflow_content = hist.GetBinContent(nbins + 1)
    overflow_error = hist.GetBinError(nbins + 1)
    last_content = hist.GetBinContent(nbins)
    last_error = hist.GetBinError(nbins)
    hist.SetBinContent(nbins, last_content + overflow_content)
    hist.SetBinError(nbins, (last_error * last_error + overflow_error * overflow_error) ** 0.5)
    hist.SetBinContent(nbins + 1, 0.0)
    hist.SetBinError(nbins + 1, 0.0)


def clip_negative_bins(hist, context, floor=0.0):
    clipped_bins = []
    for bin_idx in range(1, hist.GetNbinsX() + 1):
        content = hist.GetBinContent(bin_idx)
        if content >= floor:
            continue
        hist.SetBinContent(bin_idx, floor)
        clipped_bins.append((bin_idx, content))
    if clipped_bins:
        print(
            f"warning: clipped {len(clipped_bins)} negative bins to {floor:g} for {context}"
        )
    return clipped_bins


def load_hist(file_path, hist_name):
    if not file_path.exists():
        return None
    root_file = ROOT.TFile.Open(str(file_path), "READ")
    if not root_file or root_file.IsZombie():
        return None
    hist = root_file.Get(hist_name)
    if not hist:
        root_file.Close()
        return None
    cloned = hist.Clone(hist_name)
    cloned.SetDirectory(0)
    root_file.Close()
    merge_underflow_overflow(cloned)
    return cloned


def load_lumi_map(path):
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def load_json_map(path):
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def expand_syst_token(name, alias_map, stack=None):
    stack = [] if stack is None else list(stack)
    if name in stack:
        cycle = " -> ".join(stack + [name])
        raise RuntimeError(f"Shape systematic alias cycle detected: {cycle}")
    values = alias_map.get(name)
    if values is None:
        return [name]

    expanded = []
    seen = set()
    next_stack = stack + [name]
    for value in values:
        for item in expand_syst_token(value, alias_map, stack=next_stack):
            if item in seen:
                continue
            seen.add(item)
            expanded.append(item)
    return expanded


def resolve_shape_systs(args):
    include_shape_syst = args.include_shape_syst or args.include_btag_syst
    if args.stat_only or not include_shape_syst:
        return []

    alias_map = load_json_map(Path(args.shape_syst_alias_config).resolve())

    if args.shape_systs is not None:
        requested_tokens = list(args.shape_systs)
    elif args.shape_syst_groups is not None:
        group_map = load_json_map(Path(args.shape_syst_group_config).resolve())
        requested_tokens = []
        for group_name in args.shape_syst_groups:
            if group_name not in group_map:
                available = ", ".join(sorted(group_map))
                raise RuntimeError(
                    f"Unknown shape systematic group {group_name}. Available groups: {available}"
                )
            requested_tokens.extend(group_map[group_name])
    elif args.include_shape_syst:
        group_map = load_json_map(Path(args.shape_syst_group_config).resolve())
        if DEFAULT_FULL_SHAPE_SYST_GROUP not in group_map:
            raise RuntimeError(
                f"Default shape systematic group {DEFAULT_FULL_SHAPE_SYST_GROUP} is missing from "
                f"{Path(args.shape_syst_group_config).resolve()}"
            )
        requested_tokens = list(group_map[DEFAULT_FULL_SHAPE_SYST_GROUP])
    else:
        requested_tokens = list(args.btag_systs)

    resolved = []
    seen = set()
    for token in requested_tokens:
        for item in expand_syst_token(token, alias_map):
            if item in seen:
                continue
            seen.add(item)
            resolved.append(item)
    return resolved


PROCESS_XS_FALLBACK_MAP = {
    "ttbar": ["ttbar"],
    "TTH": ["TTHnonBB"],
    "TV": ["TZQB"],
    "TTBB": ["TTBB"],
    "TTZ": ["TTZ_low", "TTZ_high"],
    "VV": ["diBoson"],
    "VVV": ["TTVV"],
    "TTW": ["TTW"],
    "TTTT": ["TTTT"],
}


def format_xs_lnN_value(value):
    if isinstance(value, (int, float)):
        up_value = float(value)
        if up_value <= 0.0:
            raise RuntimeError(f"Invalid XS lnN value: {value}")
        return f"{up_value:g}"
    if not isinstance(value, list) or len(value) != 2:
        raise RuntimeError(
            f"XS lnN entry must be a positive number or [down, up], got: {value}"
        )
    down_value = float(value[0])
    up_value = float(value[1])
    if down_value <= 0.0 or up_value <= 0.0:
        raise RuntimeError(f"Invalid XS lnN values: {value}")
    if abs(down_value - up_value) < 1e-12:
        return f"{up_value:g}"
    return f"{down_value:g}/{up_value:g}"


def resolve_process_xs_value(process, raw_map):
    if process in raw_map:
        return format_xs_lnN_value(raw_map[process])

    fallback_keys = PROCESS_XS_FALLBACK_MAP.get(process, [])
    if not fallback_keys:
        return None

    resolved_values = []
    missing_keys = []
    for key in fallback_keys:
        if key not in raw_map:
            missing_keys.append(key)
            continue
        resolved_values.append(format_xs_lnN_value(raw_map[key]))

    if not resolved_values:
        if process == "TTBB":
            return "0.76/1.24"
        return None
    if len(set(resolved_values)) != 1:
        raise RuntimeError(
            f"Inconsistent XS lnN entries for merged process {process}: "
            + ", ".join(f"{key}={format_xs_lnN_value(raw_map[key])}" for key in fallback_keys if key in raw_map)
        )
    if missing_keys and process != "TTBB":
        raise RuntimeError(
            f"Incomplete XS lnN mapping for merged process {process}. Missing keys: {', '.join(missing_keys)}"
        )
    return resolved_values[0]


def load_xs_lnN_map(path, processes=None):
    with open(path, "r", encoding="utf-8") as handle:
        raw_map = json.load(handle)

    xs_map = {}
    for process in processes if processes is not None else PROCESSES:
        if process == "TTHH":
            continue
        resolved_value = resolve_process_xs_value(process, raw_map)
        if resolved_value is not None:
            xs_map[process] = resolved_value
    return xs_map


def sum_campaign_lumi(lumi_map, campaigns, context):
    if not campaigns:
        raise RuntimeError(f"Empty campaign list for {context}.")
    total = 0.0
    missing = []
    for campaign in campaigns:
        if campaign not in lumi_map:
            missing.append(campaign)
            continue
        total += float(lumi_map[campaign])
    if missing:
        raise RuntimeError(f"Missing lumi entries for {context}: {', '.join(missing)}")
    return total


def merge_process_hist(
    base_dir,
    category,
    process,
    campaigns,
    flat_input=False,
    hist_base_name="Signal_bestScore",
    folder_template=None,
    scale_after_merge=1.0,
    input_weight="MCTotalWeight",
):
    merged = None
    contributors = []
    hist_name = f"{process}_{hist_base_name}_{input_weight}"
    for campaign in campaigns:
        if flat_input:
            if folder_template:
                campaign_dir = base_dir / folder_template.format(category=category, campaign=campaign)
            else:
                campaign_dir = base_dir / campaign if (base_dir / campaign).is_dir() else base_dir
            file_path = campaign_dir / f"{hist_name}.root"
        else:
            file_path = category_dir(base_dir, category, campaign, folder_template=folder_template) / f"{hist_name}.root"
        hist = load_hist(file_path, hist_name)
        if hist is None:
            continue
        integral = hist.Integral()
        if integral == 0.0:
            continue
        if merged is None:
            merged = hist.Clone(hist_name)
            merged.SetDirectory(0)
        else:
            merged.Add(hist)
        contributors.append((campaign, integral))
    if merged is not None and scale_after_merge != 1.0:
        merged.Scale(scale_after_merge)
    if merged is not None:
        clip_negative_bins(merged, f"nominal {category}/{process}")
    return merged, contributors


def hist_path(base_dir, category, campaign, hist_name, flat_input=False, folder_template=None):
    if flat_input:
        if folder_template:
            campaign_dir = base_dir / folder_template.format(category=category, campaign=campaign)
        else:
            campaign_dir = base_dir / campaign if (base_dir / campaign).is_dir() else base_dir
        return campaign_dir / f"{hist_name}.root"
    return category_dir(base_dir, category, campaign, folder_template=folder_template) / f"{hist_name}.root"


def merge_process_variation_hist(
    base_dir,
    category,
    process,
    campaigns,
    syst_name,
    direction,
    flat_input=False,
    hist_base_name="Signal_bestScore",
    folder_template=None,
    input_weight="MCTotalWeight",
):
    merged = None
    used_any_variation = False
    nominal_hist_name = f"{process}_{hist_base_name}_{input_weight}"
    variation_hist_name = f"{process}_{hist_base_name}_{input_weight}_{syst_name}_{direction}"

    for campaign in campaigns:
        variation_path = hist_path(
            base_dir,
            category,
            campaign,
            variation_hist_name,
            flat_input=flat_input,
            folder_template=folder_template,
        )
        hist = load_hist(variation_path, variation_hist_name)
        if hist is not None:
            used_any_variation = True
        else:
            nominal_path = hist_path(
                base_dir,
                category,
                campaign,
                nominal_hist_name,
                flat_input=flat_input,
                folder_template=folder_template,
            )
            hist = load_hist(nominal_path, nominal_hist_name)
        if hist is None:
            continue
        if merged is None:
            merged = hist.Clone(variation_hist_name)
            merged.SetDirectory(0)
        else:
            merged.Add(hist)

    if not used_any_variation:
        return None
    if merged is not None:
        clip_negative_bins(merged, f"variation {category}/{process}/{syst_name}/{direction}")
        merged.SetName(variation_hist_name)
        merged.SetTitle(variation_hist_name)
    return merged


def ensure_can_write(path, overwrite=False, kind="output"):
    if path.exists() and not overwrite:
        raise FileExistsError(
            f"Refusing to overwrite existing {kind}: {path}. "
            "Pass --overwrite-output to allow replacing it."
        )


def write_hist(file_path, hist, overwrite=False):
    ensure_can_write(file_path, overwrite=overwrite, kind="ROOT file")
    file_path.parent.mkdir(parents=True, exist_ok=True)
    root_file = ROOT.TFile.Open(str(file_path), "RECREATE")
    root_file.cd()
    root_file.WriteObject(hist, hist.GetName())
    root_file.Close()


def build_outputs(
    base_dir,
    out_dir,
    out_tag,
    campaigns,
    categories,
    flat_input=False,
    hist_base_name="Signal_bestScore",
    shape_systs=None,
    folder_template=None,
    signal_campaigns=None,
    signal_scale_after_merge=1.0,
    overwrite_output=False,
    processes=None,
    dd_processes=None,
):
    merged_dirs = {}
    category_processes = {}
    observations = {}
    available_shape_systs = {}
    shape_systs = shape_systs or []

    processes = processes or PROCESSES
    dd_processes = set(dd_processes or [])

    for category in categories:
        if flat_input:
            merged_dir = out_dir / f"S4_dimuon_DNN_{category}_{out_tag}"
        elif category == "TTTT":
            merged_dir = out_dir / f"S4_dimuon_DNN_TTTT_region_{out_tag}"
        else:
            merged_dir = out_dir / f"S4_dimuon_DNN_{category}_region_{out_tag}"
        if merged_dir.exists() and not overwrite_output:
            raise FileExistsError(
                f"Refusing to reuse existing merged output directory: {merged_dir}. "
                "Pass --overwrite-output to allow replacing files inside it."
            )
        merged_dirs[category] = merged_dir
        category_processes[category] = []
        available_shape_systs[category] = {}
        obs = 0.0

        for process in processes:
            process_campaigns = signal_campaigns if (process == "TTHH" and signal_campaigns is not None) else campaigns
            process_scale_after_merge = signal_scale_after_merge if process == "TTHH" else 1.0
            merged_hist, contributors = merge_process_hist(
                base_dir,
                category,
                process,
                process_campaigns,
                flat_input=flat_input,
                hist_base_name=hist_base_name,
                folder_template=folder_template,
                scale_after_merge=process_scale_after_merge,
                input_weight="DDTotalWeight" if process in dd_processes else "MCTotalWeight",
            )
            if merged_hist is None:
                continue
            if merged_hist.Integral() <= 0.0:
                continue
            output_file = merged_dir / f"{process}_{hist_base_name}_MCTotalWeight.root"
            write_hist(output_file, merged_hist, overwrite=overwrite_output)
            category_processes[category].append(process)
            obs += merged_hist.Integral()
            available_shape_systs[category][process] = []

            if process == "TTHH":
                continue
            for syst_name in shape_systs:
                up_hist = merge_process_variation_hist(
                    base_dir,
                    category,
                    process,
                    campaigns,
                    syst_name,
                    "up",
                    flat_input=flat_input,
                    hist_base_name=hist_base_name,
                    folder_template=folder_template,
                    input_weight="DDTotalWeight" if process in dd_processes else "MCTotalWeight",
                )
                down_hist = merge_process_variation_hist(
                    base_dir,
                    category,
                    process,
                    campaigns,
                    syst_name,
                    "down",
                    flat_input=flat_input,
                    hist_base_name=hist_base_name,
                    folder_template=folder_template,
                    input_weight="DDTotalWeight" if process in dd_processes else "MCTotalWeight",
                )
                if up_hist is None or down_hist is None:
                    continue
                if up_hist.Integral() <= 0.0 or down_hist.Integral() <= 0.0:
                    print(
                        "warning: skip shape systematic "
                        f"{category}/{process}/{syst_name} because "
                        f"integral is non-positive after clipping "
                        f"(up={up_hist.Integral() if up_hist is not None else 'missing'}, "
                        f"down={down_hist.Integral() if down_hist is not None else 'missing'})"
                    )
                    continue
                root_file = ROOT.TFile.Open(str(output_file), "UPDATE")
                root_file.cd()
                up_hist.SetName(f"{process}_{hist_base_name}_MCTotalWeight_{syst_name}Up")
                up_hist.SetTitle(up_hist.GetName())
                down_hist.SetName(f"{process}_{hist_base_name}_MCTotalWeight_{syst_name}Down")
                down_hist.SetTitle(down_hist.GetName())
                root_file.WriteObject(up_hist, up_hist.GetName(), "Overwrite")
                root_file.WriteObject(down_hist, down_hist.GetName(), "Overwrite")
                root_file.Close()
                available_shape_systs[category][process].append(syst_name)

        observations[category] = obs

    return merged_dirs, category_processes, observations, available_shape_systs


def build_card_text(
    out_tag,
    merged_dirs,
    category_processes,
    observations,
    lumi_value,
    categories,
    hist_base_name,
    stat_only=False,
    shape_systs=None,
    available_shape_systs=None,
    xs_lnN_map=None,
    processes=None,
):
    lines = [
        f"imax {len(categories)}  number of channels",
        "jmax *  number of backgrounds",
        "kmax *  number of nuisance parameters",
        "",
    ]

    for category in categories:
        lines.append(f"shapes data_obs {category} FAKE")

    for category in categories:
        dir_path = merged_dirs[category].as_posix()
        lines.append(
            f"shapes * {category} {dir_path}/$PROCESS_{hist_base_name}_MCTotalWeight.root "
            f"$PROCESS_{hist_base_name}_MCTotalWeight "
            f"$PROCESS_{hist_base_name}_MCTotalWeight_$SYSTEMATIC"
        )

    lines.extend(
        [
            "",
            "bin           " + " ".join(categories),
            "observation   "
            + " ".join(str(observations[category]) for category in categories),
            "",
        ]
    )

    bin_tokens = []
    proc_names = []
    proc_ids = []
    rates = []
    lumi_tokens = []
    next_id = 0

    for category in categories:
        for process in category_processes[category]:
            bin_tokens.append(category)
            proc_names.append(process)
            if process == "TTHH":
                proc_ids.append("0")
            else:
                next_id += 1
                proc_ids.append(str(next_id))
            rates.append("-1")
            lumi_tokens.append(lumi_value)

    lines.append("bin           " + " ".join(bin_tokens))
    lines.append("process       " + " ".join(proc_names))
    lines.append("process       " + " ".join(proc_ids))
    lines.append("rate          " + " ".join(rates))
    lines.append("")
    if not stat_only:
        lines.append("lumi          lnN " + " ".join(lumi_tokens))
        for process_name in processes or PROCESSES:
            if process_name not in (xs_lnN_map or {}):
                continue
            token_value = xs_lnN_map[process_name]
            tokens = []
            has_effect = False
            for process in proc_names:
                if process == process_name:
                    tokens.append(token_value)
                    has_effect = True
                else:
                    tokens.append("-")
            if has_effect:
                lines.append(f"xs_{process_name:<10} lnN " + " ".join(tokens))
        for syst_name in shape_systs or []:
            tokens = []
            has_effect = False
            for category, process in zip(bin_tokens, proc_names):
                process_systs = (available_shape_systs or {}).get(category, {}).get(process, [])
                if syst_name in process_systs:
                    tokens.append("1")
                    has_effect = True
                else:
                    tokens.append("-")
            if has_effect:
                lines.append(f"{syst_name:<13} shape " + " ".join(tokens))
    lines.append("")

    return "\n".join(lines)


def main():
    args = parse_args()
    base_dir = Path(args.input_dir).resolve()
    out_dir = Path(args.output_dir).resolve()

    categories = list(args.categories)
    if "TTHH" in args.background_processes:
        raise RuntimeError("TTHH is the fixed signal and must not be listed in --background-processes.")
    if len(set(args.background_processes)) != len(args.background_processes):
        raise RuntimeError("--background-processes contains duplicate channel names.")
    processes = ["TTHH"] + list(args.background_processes)
    unknown_dd = set(args.dd_processes) - set(args.background_processes)
    if unknown_dd:
        raise RuntimeError(
            "Every --dd-processes entry must also be in --background-processes: "
            + ", ".join(sorted(unknown_dd))
        )
    if args.flat_input:
        if not args.single_region_bin:
            raise RuntimeError("--flat-input requires --single-region-bin.")
        categories = [args.single_region_bin]

    shape_systs = resolve_shape_systs(args)
    signal_campaigns = list(args.signal_campaigns) if args.signal_campaigns else None
    xs_lnN_map = {}
    if not args.stat_only and not args.disable_xs_lnN:
        xs_config_path = Path(args.xs_syst_config).resolve()
        xs_lnN_map = load_xs_lnN_map(xs_config_path, processes=processes)
    signal_scale_after_merge = 1.0
    if signal_campaigns is not None:
        lumi_map = load_lumi_map(args.lumi_config)
        target_lumi = sum_campaign_lumi(lumi_map, args.campaigns, "background/target campaigns")
        signal_lumi = sum_campaign_lumi(lumi_map, signal_campaigns, "signal campaigns")
        if signal_lumi <= 0.0:
            raise RuntimeError("Signal campaign lumi must be positive.")
        signal_scale_after_merge = target_lumi / signal_lumi

    merged_dirs, category_processes, observations, available_shape_systs = build_outputs(
        base_dir,
        out_dir,
        args.output_tag,
        args.campaigns,
        categories,
        flat_input=args.flat_input,
        hist_base_name=args.hist_base_name,
        shape_systs=shape_systs,
        folder_template=args.folder_template,
        signal_campaigns=signal_campaigns,
        signal_scale_after_merge=signal_scale_after_merge,
        overwrite_output=args.overwrite_output,
        processes=processes,
        dd_processes=args.dd_processes,
    )

    for category in categories:
        if not category_processes[category]:
            raise RuntimeError(f"Category {category} has no non-zero processes after merging.")

    card_text = build_card_text(
        args.output_tag,
        merged_dirs,
        category_processes,
        observations,
        args.lumi_lnN,
        categories,
        args.hist_base_name,
        stat_only=args.stat_only,
        shape_systs=shape_systs,
        available_shape_systs=available_shape_systs,
        xs_lnN_map=xs_lnN_map,
        processes=processes,
    )
    card_path = out_dir / f"shape_fit_{args.output_tag}.txt"
    ensure_can_write(card_path, overwrite=args.overwrite_output, kind="datacard")
    card_path.write_text(card_text, encoding="utf-8")

    print(f"merged campaigns: {', '.join(args.campaigns)}")
    if signal_campaigns is not None:
        print(
            "signal campaigns: "
            + ", ".join(signal_campaigns)
            + f" (extra scale to target lumi = {signal_scale_after_merge:.6g})"
        )
    if shape_systs:
        print(f"shape syst count: {len(shape_systs)}")
        print("shape systs: " + ", ".join(shape_systs))
    if not args.stat_only and not args.disable_xs_lnN:
        print(f"XS lnN config: {Path(args.xs_syst_config).resolve()}")
    for category in categories:
        print(f"{category}: observation={observations[category]}")
        print(f"{category}: processes={', '.join(category_processes[category])}")
    print(f"card: {card_path}")
    print(
        f"text2workspace.py {card_path.name} -o signalBestScore_shape_{args.output_tag}.root"
    )
    print(
        f"combine -M AsymptoticLimits signalBestScore_shape_{args.output_tag}.root -t -1 --expectSignal 0 -n .signalBestScoreShape_{args.output_tag}"
    )


if __name__ == "__main__":
    main()
