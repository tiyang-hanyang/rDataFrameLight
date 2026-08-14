import argparse
import json
import fnmatch
from pathlib import Path, PurePosixPath

import ROOT

'''
This script is added to jet-corrected MR, example run:
python3 ../../source/script/checkTriggerEff.py --sample-json $HOME/public/rDataFrameLight_update/source/json/samples/measurementRegion/PUJECcorrected/RunIII2024Summer24NanoAODv15_corrected.json --channels WJet*  --era RunIII2024Summer24NanoAODv15 --trigger-scheme Mu_8_27 --threshold-scheme Mar25_tuned --output-dir triggerEff_relaxed_Mar25
Could use: --selection-stage <loose/fakeable> --triggers all --variables all
'''

ROOT.EnableImplicitMT()

# assist functions for the RDF usage
ROOT.gInterpreter.Declare(
    r"""
    #include "ROOT/RVec.hxx"
    #include <cmath>

    float localDeltaPhi(float phi1, float phi2) {
        float dphi = phi1 - phi2;
        constexpr float pi = 3.14159265358979323846f;
        while (dphi > pi) dphi -= 2.0f * pi;
        while (dphi <= -pi) dphi += 2.0f * pi;
        return dphi;
    }

    float localDeltaR(float eta1, float eta2, float phi1, float phi2) {
        const float deta = eta1 - eta2;
        const float dphi = localDeltaPhi(phi1, phi2);
        return std::sqrt(deta * deta + dphi * dphi);
    }

    ROOT::VecOps::RVec<float> allJetDRFromMuon(ROOT::VecOps::RVec<float> Jet_eta, float Muon_eta, ROOT::VecOps::RVec<float> Jet_phi, float Muon_phi) {
        ROOT::VecOps::RVec<float> DR;
        int nJet = Jet_eta.size();
        for (int i = 0; i < nJet; i++) {
            DR.push_back(localDeltaR(Jet_eta[i], Muon_eta, Jet_phi[i], Muon_phi));
        }
        return DR;
    }

    ROOT::VecOps::RVec<float> getMuonMatchJetDR(const ROOT::VecOps::RVec<short>& muon_jetIdx, const ROOT::VecOps::RVec<float>& jet_eta, const ROOT::VecOps::RVec<float>& jet_phi, const ROOT::VecOps::RVec<float>& muon_eta, const ROOT::VecOps::RVec<float>& muon_phi) {
        ROOT::VecOps::RVec<float> muonjetdr;
        int nMuon = muon_jetIdx.size();
        for (int i = 0; i < nMuon; i++) {
            if (muon_jetIdx[i] == -1) {
                muonjetdr.push_back(999.f);
                continue;
            }
            muonjetdr.push_back(localDeltaR(muon_eta[i], jet_eta[muon_jetIdx[i]], muon_phi[i], jet_phi[muon_jetIdx[i]]));
        }
        return muonjetdr;
    }

    ROOT::VecOps::RVec<float> getMuonMatchJetPt(const ROOT::VecOps::RVec<short>& muon_jetIdx, const ROOT::VecOps::RVec<float>& jet_pt) {
        ROOT::VecOps::RVec<float> muonjetpt;
        int nMuon = muon_jetIdx.size();
        for (int i = 0; i < nMuon; i++) {
            if (muon_jetIdx[i] == -1) {
                muonjetpt.push_back(-999.f);
                continue;
            }
            muonjetpt.push_back(jet_pt[muon_jetIdx[i]]);
        }
        return muonjetpt;
    }

    ROOT::VecOps::RVec<short> ensureMuonJetGood(const ROOT::VecOps::RVec<short>& muon_jetIdx, const ROOT::VecOps::RVec<bool>& GoodJetCond) {
        ROOT::VecOps::RVec<short> Muon_jetIdxGood;
        int nMuon = muon_jetIdx.size();
        for (int i = 0; i < nMuon; i++) {
            if (muon_jetIdx[i] == -1) {
                Muon_jetIdxGood.push_back(-1);
                continue;
            }
            if (GoodJetCond[muon_jetIdx[i]]) {
                Muon_jetIdxGood.push_back(muon_jetIdx[i]);
            } else {
                Muon_jetIdxGood.push_back(-1);
            }
        }
        return Muon_jetIdxGood;
    }

    """
)

# Trigger thresholds are expressed only once. The pass expression is derived
# from the trigger name itself, while the denominator reuses the same threshold.
TRIGGER_THRESHOLDS = {
    "Run2_prescaled": {
        "HLT_Mu3_PFJet40": "(leadingMuonRecoPt>3 && leadingMuonConePt<32 && leadingIsoJetPt>45)",
        "HLT_Mu8": "(leadingMuonRecoPt>8 && leadingMuonConePt<100)",
        "HLT_Mu17": "(leadingMuonRecoPt>17 && leadingMuonConePt<100 && leadingMuonConePt>32)",
        "HLT_Mu20": "(leadingMuonRecoPt>20 && leadingMuonConePt<100 && leadingMuonConePt>32)",
        "HLT_Mu27": "(leadingMuonRecoPt>27 && leadingMuonConePt<100 && leadingMuonConePt>45)",
    },
    "Mar25_tuned": {
        "HLT_Mu3_PFJet40": "(leadingMuonRecoPt>3 && leadingMuonConePt<32 && leadingIsoJetPt>45)",
        "HLT_Mu8": "(leadingMuonRecoPt>10 && leadingMuonConePt<100)",
        "HLT_Mu17": "(leadingMuonRecoPt>17 && leadingMuonConePt<100 && leadingMuonConePt>32)",
        "HLT_Mu20": "(leadingMuonRecoPt>30 && leadingMuonConePt<100 && leadingMuonConePt>40)",
        "HLT_Mu27": "(leadingMuonRecoPt>40 && leadingMuonConePt<100 && leadingMuonConePt>55)",
    },
}

TRIGGER_SCHEMES = {
    "Run2_prescaled_trig": ["HLT_Mu3_PFJet40", "HLT_Mu8", "HLT_Mu17", "HLT_Mu20", "HLT_Mu27"],
    "Mu_8_27": ["HLT_Mu8", "HLT_Mu17", "HLT_Mu20", "HLT_Mu27"],
}

PRESCALE_VALUES = {
    "HLT_Mu3_PFJet40": 22000.0,
    "HLT_Mu8": 8800.0,
    "HLT_Mu17": 325.0,
    "HLT_Mu20": 1050.0,
    "HLT_Mu27": 500.0,
}

VARIABLE_CONFIGS = {
    "leadingIsoJetPt": {
        "title": "leading isolated jet p_{T} [GeV]",
        "nbins": 30,
        "xmin": 0.0,
        "xmax": 300.0,
    },
    "leadingMuonConePt": {
        "title": "leading muon cone p_{T} [GeV]",
        "nbins": 20,
        "xmin": 0.0,
        "xmax": 100.0,
    },
    "leadingMuonPt": {
        "title": "leading muon reco p_{T} [GeV]",
        "nbins": 20,
        "xmin": 0.0,
        "xmax": 100.0,
    },
}


def parse_args():
    parser = argparse.ArgumentParser(description="Check MR trigger efficiency with configurable muon selection and threshold schemes.")
    parser.add_argument("--sample-json", help="Sample json path used in MR mode.")
    parser.add_argument("--channels", nargs="+", help="One or more channels or wildcard patterns to combine.")
    parser.add_argument(
        "--event-fraction-per-channel",
        type=float,
        default=1.0,
        help="Fraction of events to read from each merged channel file. Must satisfy 0 < fraction <= 1. Default uses all events.",
    )
    parser.add_argument("--era", help="Era key for lumi lookup, e.g. RunIII2024Summer24NanoAODv15.")
    parser.add_argument("--output-dir", default="triggerEff_MR")
    parser.add_argument(
        "--triggers",
        nargs="+",
        default=["all"],
        help="Trigger list. Use explicit names, 'single', 'or', or 'all'. The trigger-scheme name itself is used for the full OR.",
    )
    parser.add_argument(
        "--threshold-scheme",
        default="Mar25_tuned",
        choices=sorted(TRIGGER_THRESHOLDS.keys()),
        help="Threshold scheme used for trigger denominator/pass definitions.",
    )
    parser.add_argument(
        "--trigger-scheme",
        default="Mu_8_27",
        choices=sorted(TRIGGER_SCHEMES.keys()),
        help="Trigger menu scheme used to define which single and combined trigger curves are produced.",
    )
    parser.add_argument(
        "--selection-stage",
        default="fakeable_muon",
        choices=["loose_muon", "fakeable_muon", "all"],
        help="Object selection used before measuring trigger efficiency.",
    )
    parser.add_argument(
        "--variables",
        nargs="+",
        default=["all"],
        help="Variable list. Use explicit names, or 'all' for every supported variable.",
    )
    parser.add_argument("--threads", type=int, default=0, help="If >0, call ROOT.EnableImplicitMT(threads).")
    parser.add_argument("--verbose", action="store_true", help="Print resolved options, input files, and processing progress.")
    parser.add_argument(
        "--notOverWriting",
        action="store_true",
        help="Do not overwrite existing output plots. Skip them instead.",
    )
    return parser.parse_args()


def load_json(path):
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def log_verbose(enabled, message):
    if enabled:
        print(f"[checkTriggerEff] {message}")


def join_sample_path(base_dir, file_name):
    if "/" in base_dir:
        return str(PurePosixPath(base_dir) / file_name)
    return str(Path(base_dir) / file_name)


def load_xs_lumi():
    source_root = Path(__file__).resolve().parents[1]
    xs_map = load_json(source_root / "json" / "XS" / "Run3.json")
    lumi_map = load_json(source_root / "json" / "Lumi" / "Run3.json")
    return xs_map, lumi_map


def build_channel_file_map(sample_json_path, channels):
    payload = load_json(sample_json_path)
    dir_map = payload.get("dir", {})
    file_map = payload.get("file", {})
    available_channels = sorted(file_map.keys())
    expanded_channels = []
    for channel_pattern in channels:
        matched = [name for name in available_channels if fnmatch.fnmatch(name, channel_pattern)]
        if not matched:
            raise RuntimeError(f"Channel pattern {channel_pattern} did not match anything in {sample_json_path}")
        expanded_channels.extend(matched)
    expanded_channels = list(dict.fromkeys(expanded_channels))
    records = {}
    for channel in expanded_channels:
        if channel not in dir_map or channel not in file_map:
            raise RuntimeError(f"Channel {channel} not found in sample json: {sample_json_path}")
        files = [join_sample_path(dir_map[channel], file_name) for file_name in file_map[channel]]
        if not files:
            raise RuntimeError(f"Channel {channel} has no files in sample json.")
        if len(files) != 1:
            raise RuntimeError(f"Channel {channel} is expected to have exactly one merged input file, found {len(files)}")
        records[channel] = files[0]
    return records


def get_fakeable_config(threshold_scheme):
    if threshold_scheme == "Mar25_tuned":
        return "(Muon_promptMVA<0.64 && (Muon_jetDF_corr<0.2480) && Muon_jetRelIso < 0.5)"
    return "(Muon_promptMVA<0.64 && Muon_jetDF_failMVA_corr && Muon_jetRelIso < 0.5)"


def resolve_trigger_requests(trigger_requests, trigger_scheme):
    single_triggers = TRIGGER_SCHEMES[trigger_scheme]
    aggregate_name = trigger_scheme
    if trigger_requests == ["all"]:
        return single_triggers + [aggregate_name]
    resolved = []
    for item in trigger_requests:
        if item == "all":
            resolved.extend(single_triggers)
            resolved.append(aggregate_name)
        elif item == "single":
            resolved.extend(single_triggers)
        elif item == "or":
            resolved.append(aggregate_name)
        else:
            resolved.append(item)
    return list(dict.fromkeys(resolved))


def build_trigger_configs(requests, trigger_scheme, threshold_scheme):
    def build_single_trigger_config(trigger_name):
        if trigger_name not in TRIGGER_SCHEMES[trigger_scheme]:
            raise RuntimeError(f"Unsupported trigger {trigger_name} for trigger scheme {trigger_scheme}")
        if trigger_name not in TRIGGER_THRESHOLDS[threshold_scheme]:
            raise RuntimeError(f"Trigger {trigger_name} is missing from threshold scheme {threshold_scheme}")
        threshold_expr = TRIGGER_THRESHOLDS[threshold_scheme][trigger_name]
        return {
            "name": trigger_name,
            "mode": "boolean",
            "denom": threshold_expr,
            "pass": f"({trigger_name} && {threshold_expr})",
        }

    def build_trigger_or_config(trigger_names, name):
        missing = [trigger for trigger in trigger_names if trigger not in TRIGGER_THRESHOLDS[threshold_scheme]]
        if missing:
            raise RuntimeError(f"Triggers {missing} are missing from threshold scheme {threshold_scheme}")
        denom_terms = [f"({TRIGGER_THRESHOLDS[threshold_scheme][trigger]})" for trigger in trigger_names]
        pass_terms = [f"pass_{trigger}" for trigger in trigger_names]
        denom_flag_terms = [f"denom_{trigger}" for trigger in trigger_names]
        weight_terms = [f"(1 - {pass_name}/{PRESCALE_VALUES[trigger]})" for pass_name, trigger in zip(pass_terms, trigger_names)]
        denom_weight_terms = [f"(1 - {denom_name}/{PRESCALE_VALUES[trigger]})" for denom_name, trigger in zip(denom_flag_terms, trigger_names)]
        return {
            "name": name,
            "mode": "weighted_or",
            "denom": "(" + " || ".join(denom_terms) + ")",
            "denom_weight_expr": "1.0 - " + " * ".join(denom_weight_terms),
            "weight_expr": "1.0 - " + " * ".join(weight_terms),
        }

    configs = []
    for trigger_name in resolve_trigger_requests(requests, trigger_scheme):
        if trigger_name == trigger_scheme:
            configs.append(build_trigger_or_config(TRIGGER_SCHEMES[trigger_scheme], trigger_name))
        else:
            configs.append(build_single_trigger_config(trigger_name))
    return configs


def resolve_variable_requests(requests):
    resolved = list(VARIABLE_CONFIGS.keys()) if "all" in requests else list(dict.fromkeys(requests))
    for item in resolved:
        if item not in VARIABLE_CONFIGS:
            raise RuntimeError(f"Unsupported variable: {item}")
    return resolved


def sanitize_label(text):
    return text.replace("*", "_all_").replace("?", "_q_")


def compute_genweight_sum(file_path):
    root_file = ROOT.TFile.Open(file_path, "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Failed to open input file: {file_path}")
    hist = root_file.Get("genWeightSum")
    if hist is None:
        root_file.Close()
        raise RuntimeError(f"Missing genWeightSum histogram in file: {file_path}")
    total = float(hist.Integral())
    root_file.Close()
    return total


def resolve_selection_stages(selection_stage):
    if selection_stage == "all":
        return ["loose_muon", "fakeable_muon"]
    return [selection_stage]


def define_common_object_columns(rdf):
    # Shared columns after the corrected-jet step. Stage-specific muon counting and
    # recoil-jet selection are applied later on branch-specific RDF nodes.
    rdf = rdf.Filter("JVMweight>0")
    rdf = rdf.Define("GoodJetForMuonMatch", "JetIdTight && (Jet_pt_JEC>15.0)")
    rdf = rdf.Define("Muon_jetIdxGood", "ensureMuonJetGood(Muon_jetIdx, GoodJetForMuonMatch)")
    rdf = rdf.Define("Muon_jetPt", "getMuonMatchJetPt(Muon_jetIdxGood, Jet_pt_JEC)")
    rdf = rdf.Define("Muon_jetDR", "getMuonMatchJetDR(Muon_jetIdxGood, Jet_eta, Jet_phi, Muon_eta, Muon_phi)")
    rdf = rdf.Define("Muon_conePt_fail", "(Muon_jetDR<=0.4)*0.9*Muon_jetPt + (Muon_jetDR>0.4)*Muon_pt/(1+Muon_miniPFRelIso_all)")
    rdf = rdf.Define("Muon_jetRelIso_corr", "(Muon_jetDR<0.4) * (Muon_jetPt /Muon_pt - 1.0) + (Muon_jetDR>=0.4) * Muon_pfRelIso04_all")
    rdf = rdf.Define("Muon_jetDF_corr", "(Muon_jetDR<0.4) * Muon_jetDF")
    rdf = rdf.Define("passMVA", "(Muon_promptMVA>=0.64 && Muon_jetDF_corr<0.2480)")
    rdf = rdf.Define("Muon_conePt", "passMVA*Muon_pt + (!passMVA)*Muon_conePt_fail")
    rdf = rdf.Define(
        "isLooseMuon",
        "(Muon_conePt>15.0) && (abs(Muon_eta)<2.4) && (abs(Muon_dxy)<0.05) && (abs(Muon_dz)<0.1) && (abs(Muon_sip3d)<8) && Muon_mediumId && (Muon_miniPFRelIso_all<0.4)"
    )
    return rdf


def define_fakeable_muon_columns(rdf, threshold_scheme):
    # Full fakeable selection used in the final MR stage. Loose already carries
    # the eta, IP, mediumId, isolation, and conePt requirements.
    rdf = rdf.Define("Muon_jetDF_threshold_high_corr", "(Muon_conePt>=45.0) && (Muon_jetDF_corr<0.0485)")
    rdf = rdf.Define(
        "Muon_jetDF_threshold_medium_corr",
        "(Muon_conePt<45.0) && (Muon_conePt>20.0) && (Muon_jetDF_corr<(45-Muon_conePt)*0.00798+0.0485)",
    )
    rdf = rdf.Define("Muon_jetDF_threshold_low_corr", "(Muon_conePt<=20.0) && (Muon_jetDF_corr<0.2480)")
    rdf = rdf.Define(
        "Muon_jetDF_failMVA_corr",
        "(Muon_jetDF_threshold_high_corr || Muon_jetDF_threshold_medium_corr || Muon_jetDF_threshold_low_corr)",
    )
    rdf = rdf.Define("failMVA", get_fakeable_config(threshold_scheme))
    rdf = rdf.Define("isFakeableMuon", "isLooseMuon && (passMVA || failMVA)")
    return rdf


def define_event_level_columns(rdf):
    # Event-level quantities are derived after selecting the stage-specific muon set.
    rdf = rdf.Define("nSelectedMuon", "Nonzero(selectedMuon).size()").Filter("nSelectedMuon == 1")
    rdf = rdf.Define("leadingMuonIdx", "Nonzero(selectedMuon)[0]")
    rdf = rdf.Define("leadingMuonPhi", "Muon_phi[leadingMuonIdx]")
    rdf = rdf.Define("leadingMuonEta", "Muon_eta[leadingMuonIdx]")
    rdf = rdf.Define("leadingMuonRecoPt", "Muon_pt[leadingMuonIdx]")
    rdf = rdf.Define("leadingMuonConePt", "Muon_conePt[leadingMuonIdx]")
    rdf = rdf.Define("leadingMuonPt", "Muon_pt[leadingMuonIdx]")
    rdf = rdf.Define("Jet_drFromLeadingMuon", "allJetDRFromMuon(Jet_eta, leadingMuonEta, Jet_phi, leadingMuonPhi)")
    rdf = rdf.Define("isIsoJet", "JetIdTight && (Jet_pt_JEC>30.0) && (Jet_drFromLeadingMuon>0.7)")
    rdf = rdf.Define("nIsoJet", "Nonzero(isIsoJet).size()").Filter("nIsoJet>0")
    rdf = rdf.Define("leadingIsoJetIdx", "Nonzero(isIsoJet)[0]")
    rdf = rdf.Define("leadingIsoJetPt", "Jet_pt_JEC[leadingIsoJetIdx]")
    return rdf


def build_stage_rdf(common_rdf, selection_stage, threshold_scheme):
    if selection_stage == "fakeable_muon":
        rdf = define_fakeable_muon_columns(common_rdf, threshold_scheme)
        rdf = rdf.Define("selectedMuon", "isFakeableMuon")
    else:
        rdf = common_rdf.Define("selectedMuon", "isLooseMuon")
    return define_event_level_columns(rdf)


def hist_model(name, cfg):
    return ROOT.RDF.TH1DModel(name, "", cfg["nbins"], cfg["xmin"], cfg["xmax"])


def merge_overflow(hist):
    nbins = hist.GetNbinsX()
    last = nbins
    hist.SetBinContent(last, hist.GetBinContent(last) + hist.GetBinContent(nbins + 1))
    hist.SetBinError(last, (hist.GetBinError(last) ** 2 + hist.GetBinError(nbins + 1) ** 2) ** 0.5)
    hist.SetBinContent(nbins + 1, 0.0)
    hist.SetBinError(nbins + 1, 0.0)
    return hist


def make_efficiency_hist(channel_hists, out_name):
    num_total = None
    den_total = None
    for num_hist, den_hist in channel_hists:
        if num_total is None:
            num_total = num_hist.Clone(out_name + "_num")
            den_total = den_hist.Clone(out_name + "_den")
            num_total.SetDirectory(0)
            den_total.SetDirectory(0)
        else:
            num_total.Add(num_hist)
            den_total.Add(den_hist)
    eff = num_total.Clone(out_name)
    eff.SetDirectory(0)
    eff.Divide(num_total, den_total, 1.0, 1.0, "B")
    return eff


def draw_efficiency(hist, title_lines, x_title, output_path):
    hist.SetStats(0)
    hist.SetTitle("")
    hist.GetYaxis().SetTitle("Trigger eff")
    hist.GetXaxis().SetTitle(x_title)
    hist.SetLineWidth(2)
    hist.SetMinimum(0.0)
    hist.SetMaximum(1.2)

    canvas = ROOT.TCanvas("c1", "c1", 700, 700)
    canvas.SetLeftMargin(0.14)
    hist.Draw("HIST E")

    latex = ROOT.TLatex()
    latex.SetTextSize(0.03)
    latex.SetTextFont(42)
    latex.SetTextAlign(12)
    latex.SetNDC(True)

    y = 0.88
    for line in title_lines:
        latex.DrawLatex(0.18, y, line)
        y -= 0.045

    canvas.SaveAs(str(output_path))


def run_mr_mode(args):
    if not args.sample_json or not args.channels or not args.era:
        raise RuntimeError("This script requires --sample-json, --channels, and --era.")
    if args.event_fraction_per_channel <= 0.0 or args.event_fraction_per_channel > 1.0:
        raise RuntimeError("--event-fraction-per-channel must satisfy 0 < fraction <= 1.")

    log_verbose(args.verbose, f"sample json: {args.sample_json}")
    log_verbose(args.verbose, f"channels request: {args.channels}")
    log_verbose(args.verbose, f"era: {args.era}")
    log_verbose(args.verbose, f"trigger scheme: {args.trigger_scheme}")
    log_verbose(args.verbose, f"threshold scheme: {args.threshold_scheme}")
    log_verbose(args.verbose, f"selection stage: {args.selection_stage}")
    log_verbose(args.verbose, f"variables request: {args.variables}")
    log_verbose(args.verbose, f"event fraction per channel: {args.event_fraction_per_channel}")
    log_verbose(args.verbose, f"output dir: {args.output_dir}")

    xs_map, lumi_map = load_xs_lumi()
    if args.era not in lumi_map:
        raise RuntimeError(f"Era {args.era} not found in json/Lumi/Run3.json")
    lumi = float(lumi_map[args.era])
    channel_files = build_channel_file_map(args.sample_json, args.channels)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    channel_label = "_".join(sanitize_label(item) for item in args.channels)
    trigger_cfgs = build_trigger_configs(args.triggers, args.trigger_scheme, args.threshold_scheme)
    variable_names = resolve_variable_requests(args.variables)
    selection_stages = resolve_selection_stages(args.selection_stage)
    log_verbose(args.verbose, f"resolved selection stages: {selection_stages}")
    log_verbose(args.verbose, f"resolved triggers: {[cfg['name'] for cfg in trigger_cfgs]}")
    log_verbose(args.verbose, f"resolved variables: {variable_names}")
    hist_pairs = {}
    for stage in selection_stages:
        for trigger_cfg in trigger_cfgs:
            for variable_name in variable_names:
                hist_pairs[(stage, trigger_cfg["name"], variable_name)] = []

    threshold_map = TRIGGER_THRESHOLDS[args.threshold_scheme]
    for channel, file_path in channel_files.items():
        log_verbose(args.verbose, f"processing channel {channel}")
        log_verbose(args.verbose, f"input file: {file_path}")
        if channel not in xs_map:
            raise RuntimeError(f"Channel {channel} not found in json/XS/Run3.json")
        gen_weight_sum = compute_genweight_sum(file_path)
        if gen_weight_sum == 0.0:
            raise RuntimeError(f"genWeight sum is zero for channel {channel}")
        scale = float(xs_map[channel]) * lumi * 1000.0 / gen_weight_sum
        log_verbose(
            args.verbose,
            f"channel {channel}: XS={float(xs_map[channel]):.6g}, lumi={lumi:.6g}, genWeightSum={gen_weight_sum:.6g}, scale={scale:.6g}",
        )

        common_rdf = ROOT.RDataFrame("Events", file_path)
        if args.event_fraction_per_channel < 1.0:
            n_events = int(common_rdf.Count().GetValue())
            n_keep = max(1, int(n_events * args.event_fraction_per_channel))
            log_verbose(args.verbose, f"channel {channel}: total events={n_events}, keeping first {n_keep}")
            common_rdf = common_rdf.Range(n_keep)
        else:
            log_verbose(args.verbose, f"channel {channel}: using all events")
        log_verbose(args.verbose, f"channel {channel}: building common object columns")
        common_rdf = define_common_object_columns(common_rdf)
        stage_rdfs = {}
        for stage in selection_stages:
            log_verbose(args.verbose, f"channel {channel}: building stage rdf for {stage}")
            stage_rdf = build_stage_rdf(common_rdf, stage, args.threshold_scheme)
            stage_rdf = stage_rdf.Define("evtWeight", f"({scale}) * genWeight")
            for single_trigger in TRIGGER_SCHEMES[args.trigger_scheme]:
                threshold_expr = threshold_map[single_trigger]
                stage_rdf = stage_rdf.Define(f"denom_{single_trigger}", threshold_expr)
                stage_rdf = stage_rdf.Define(
                    f"pass_{single_trigger}",
                    f"({single_trigger} && {threshold_expr})",
                )
            stage_rdfs[stage] = stage_rdf

        for stage, stage_rdf in stage_rdfs.items():
            for trigger_cfg in trigger_cfgs:
                log_verbose(args.verbose, f"channel {channel}: stage={stage}, trigger={trigger_cfg['name']}")
                denom_rdf = stage_rdf.Filter(trigger_cfg["denom"])
                if trigger_cfg["mode"] == "weighted_or":
                    denom_rdf = denom_rdf.Define("triggerDenominatorWeight", trigger_cfg["denom_weight_expr"])
                    denom_rdf = denom_rdf.Define("triggerDecisionWeight", trigger_cfg["weight_expr"])
                    denom_rdf = denom_rdf.Define("evtWeightDenomTrigger", "evtWeight * triggerDenominatorWeight")
                    denom_rdf = denom_rdf.Define("evtWeightTrigger", "evtWeight * triggerDecisionWeight")
                for variable_name in variable_names:
                    log_verbose(args.verbose, f"channel {channel}: stage={stage}, trigger={trigger_cfg['name']}, variable={variable_name}")
                    var_cfg = VARIABLE_CONFIGS[variable_name]
                    if trigger_cfg["mode"] == "weighted_or":
                        num_hist = denom_rdf.Histo1D(
                            hist_model(f"h_num_{trigger_cfg['name']}_{variable_name}_{channel}_{stage}", var_cfg),
                            variable_name,
                            "evtWeightTrigger",
                        ).GetValue()
                        den_hist = denom_rdf.Histo1D(
                            hist_model(f"h_den_{trigger_cfg['name']}_{variable_name}_{channel}_{stage}", var_cfg),
                            variable_name,
                            "evtWeightDenomTrigger",
                        ).GetValue()
                    else:
                        num_hist = denom_rdf.Filter(trigger_cfg["pass"]).Histo1D(
                            hist_model(f"h_num_{trigger_cfg['name']}_{variable_name}_{channel}_{stage}", var_cfg),
                            variable_name,
                            "evtWeight",
                        ).GetValue()
                        den_hist = denom_rdf.Histo1D(
                            hist_model(f"h_den_{trigger_cfg['name']}_{variable_name}_{channel}_{stage}", var_cfg),
                            variable_name,
                            "evtWeight",
                        ).GetValue()
                    num_hist.SetDirectory(0)
                    den_hist.SetDirectory(0)
                    merge_overflow(num_hist)
                    merge_overflow(den_hist)
                    hist_pairs[(stage, trigger_cfg["name"], variable_name)].append((num_hist, den_hist))

    for stage in selection_stages:
        stage_label = stage.replace("_", "-")
        for trigger_cfg in trigger_cfgs:
            for variable_name in variable_names:
                var_cfg = VARIABLE_CONFIGS[variable_name]
                output_path = output_dir / f"{channel_label}_{stage_label}_{trigger_cfg['name']}_{variable_name}.png"
                if args.notOverWriting and output_path.exists():
                    print(f"[checkTriggerEff] skip existing plot: {output_path}")
                    continue
                log_verbose(args.verbose, f"writing plot: {output_path}")
                eff_hist = make_efficiency_hist(
                    hist_pairs[(stage, trigger_cfg["name"], variable_name)],
                    f"eff_{trigger_cfg['name']}_{variable_name}_{channel_label}_{stage}",
                )
                draw_efficiency(
                    eff_hist,
                    [
                        f"Channels: {', '.join(args.channels)}",
                        f"Selection: {stage.replace('_', ' ')}",
                        f"Trigger: {trigger_cfg['name']}",
                        f"Trigger scheme: {args.trigger_scheme}",
                        f"Threshold scheme: {args.threshold_scheme}",
                    ],
                    var_cfg["title"],
                    output_path,
                )


def main():
    args = parse_args()
    if args.threads and args.threads > 0:
        ROOT.EnableImplicitMT(args.threads)

    run_mr_mode(args)


if __name__ == "__main__":
    main()
