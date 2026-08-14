import argparse
import json
from pathlib import Path

import numpy as np
import uproot

PAIR_INDEX = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
REQUIRED_BRANCHES = [
    "run",
    "luminosityBlock",
    "event",
    "Jet_pt_JEC",
    "Jet_eta",
    "Jet_phi",
    "Jet_mass_JEC",
    "leadingMuonIdx",
    "subleadingMuonIdx",
    "Muon_charge",
    "Muon_eta",
    "Muon_phi",
    "Muon_mass",
    "Muon_miniPFRelIso_all",
    "Muon_jetRelIso",
    "Muon_jetDF",
]
OPTIONAL_JET_SCORE_BRANCHES = ["Jet_btagUParTAK4B", "Jet_btagPNetB"]
OPTIONAL_MUON_PROMPT_BRANCHES = ["Muon_promptMVA", "Muon_mvaTTH"]
OPTIONAL_PT_BRANCHES = ["Muon_pt_Rcorr", "Muon_pt", "Muon_pt_Rscale"]
OPTIONAL_WEIGHT_BRANCHES = ["genWeight", "PUWeight"]
JET_SELECTION_BRANCHES = ["Jet_rawFactor", "JetIdTight", "Jet_drFromMuon"]
JET_BRANCH_ALIASES = {
    "Jet_pt_JEC": "Jet_pt",
    "Jet_mass_JEC": "Jet_mass",
}
BJET_NODE_FIELD_NAMES = []
for jet_index in range(1, 5):
    BJET_NODE_FIELD_NAMES.extend(
        [
            f"jet{jet_index}_px",
            f"jet{jet_index}_py",
            f"jet{jet_index}_pz",
            f"jet{jet_index}_mass",
            f"jet{jet_index}_btag",
            f"jet{jet_index}_dr_mu1",
            f"jet{jet_index}_dr_mu2",
        ]
    )

FEATURE_FIELD_NAMES = (
    BJET_NODE_FIELD_NAMES
    + [f"bjet_dr_{pair[0]}{pair[1]}" for pair in PAIR_INDEX]
    + [f"bjet_m_{pair[0]}{pair[1]}" for pair in PAIR_INDEX]
    + [
        "lead_mu_px",
        "lead_mu_py",
        "lead_mu_pz",
        "sublead_mu_px",
        "sublead_mu_py",
        "sublead_mu_pz",
        "lead_mu_miniPFRelIso_all",
        "sublead_mu_miniPFRelIso_all",
        "lead_mu_jetRelIso",
        "sublead_mu_jetRelIso",
        "lead_mu_jetDF",
        "sublead_mu_jetDF",
        "lead_mu_promptMVA",
        "sublead_mu_promptMVA",
        "dimuon_dr",
        "dimuon_mass",
        "nGoodJet",
        "nBJet",
        "sumjet_pt",
        "sumbjet_pt",
        "jet_cent",
    ]
)
ALL_FIELD_NAMES = FEATURE_FIELD_NAMES + ["weight", "isSameSign"]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract simple-DNN signal-classification inputs from NanoAOD-like ROOT files."
    )
    parser.add_argument("--input", required=True, help="Input ROOT file.")
    parser.add_argument("--output", required=True, help="Output .npy path.")
    parser.add_argument("--sample-name", default=None, help="Sample name used for label inference.")
    parser.add_argument("--era", required=True, help="Campaign key used to resolve luminosity.")
    parser.add_argument("--tree-name", default="Events", help="Tree name.")
    parser.add_argument("--verbose", type=int, default=1, help="Set to 0 to disable normalization debug prints.")
    parser.add_argument("--jet-pt-shift-branch", default=None, help="JES shift branch, e.g. CMS_scale_j_FlavorQCD.")
    parser.add_argument("--jet-pt-shift-direction", choices=["up", "down"], default="up")
    parser.add_argument("--jer-direction", choices=["up", "down"], default=None, help="Use Jet_JER_corr_{up/down} to vary Jet_pt.")
    parser.add_argument(
        "--muon-pt-branch",
        default=None,
        help=(
            "Optional override for the muon pt branch used in muon four-vectors and dimuon mass, "
            "e.g. Muon_pt_Rcorr, Muon_pt_Rscale_up, Muon_pt_Rcorr_resolup."
        ),
    )
    parser.add_argument("--jet-pt-threshold", type=float, default=30.0)
    parser.add_argument("--btag-threshold", type=float, default=None)
    parser.add_argument("--scale-jet-mass-with-pt", action="store_true")
    return parser.parse_args()


def resolve_sample_name(sample_name, input_path, xs_map):
    if sample_name is not None:
        return sample_name

    file_stem = Path(input_path).stem
    matched_keys = [key for key in xs_map if key in file_stem]
    if not matched_keys:
        raise RuntimeError(
            f"Unable to resolve sample name from file name: {file_stem}. "
            "Pass --sample-name explicitly."
        )
    matched_keys.sort(key=len, reverse=True)
    return matched_keys[0]


def delta_phi(phi1, phi2):
    dphi = phi1 - phi2
    return (dphi + np.pi) % (2 * np.pi) - np.pi


def delta_r(eta1, phi1, eta2, phi2):
    deta = eta1 - eta2
    dphi = delta_phi(phi1, phi2)
    return np.sqrt(deta * deta + dphi * dphi)


def energy_from_pt_eta_mass(pt, eta, mass):
    return np.sqrt((pt * np.cosh(eta)) ** 2 + mass**2)


def four_vector_from_pt_eta_phi_mass(pt, eta, phi, mass):
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    energy = energy_from_pt_eta_mass(pt, eta, mass)
    return px, py, pz, energy


def invariant_mass(pt1, eta1, phi1, mass1, pt2, eta2, phi2, mass2):
    px1, py1, pz1, energy1 = four_vector_from_pt_eta_phi_mass(pt1, eta1, phi1, mass1)
    px2, py2, pz2, energy2 = four_vector_from_pt_eta_phi_mass(pt2, eta2, phi2, mass2)
    m2 = (energy1 + energy2) ** 2 - (
        (px1 + px2) ** 2 + (py1 + py2) ** 2 + (pz1 + pz2) ** 2
    )
    return np.sqrt(np.maximum(m2, 0.0))


def load_arrays(
    file_name,
    tree_name,
    jet_pt_shift_branch=None,
    jer_direction=None,
    muon_pt_branch=None,
):
    with uproot.open(file_name) as root_file:
        tree = root_file[tree_name]
        all_branches = {branch.name for branch in tree.branches}

        missing = [name for name in REQUIRED_BRANCHES if name not in all_branches]
        if missing:
            raise RuntimeError("Missing required branches: " + ", ".join(missing))

        jet_score_branch = next((name for name in OPTIONAL_JET_SCORE_BRANCHES if name in all_branches), None)
        if jet_score_branch is None:
            raise RuntimeError(
                "Unable to find a jet score branch from: " + ", ".join(OPTIONAL_JET_SCORE_BRANCHES)
            )

        muon_prompt_branch = next((name for name in OPTIONAL_MUON_PROMPT_BRANCHES if name in all_branches), None)
        if muon_prompt_branch is None:
            raise RuntimeError(
                "Unable to find a muon prompt-MVA branch from: " + ", ".join(OPTIONAL_MUON_PROMPT_BRANCHES)
            )

        if muon_pt_branch is not None:
            if muon_pt_branch not in all_branches:
                raise RuntimeError(f"Missing requested muon pt branch: {muon_pt_branch}")
            pt_branch = muon_pt_branch
        else:
            pt_branch = next((name for name in OPTIONAL_PT_BRANCHES if name in all_branches), None)
            if pt_branch is None:
                raise RuntimeError(
                    "Unable to find a muon pt branch from: " + ", ".join(OPTIONAL_PT_BRANCHES)
                )

        branches = REQUIRED_BRANCHES + [jet_score_branch, muon_prompt_branch, pt_branch]
        if jet_pt_shift_branch is not None:
            if jet_pt_shift_branch not in all_branches:
                raise RuntimeError(f"Missing requested jet pt shift branch: {jet_pt_shift_branch}")
            branches.append(jet_pt_shift_branch)
        if jer_direction is not None:
            jer_branches = ["Jet_JER_corr", f"Jet_JER_corr_{jer_direction}"]
            missing_jer = [name for name in jer_branches if name not in all_branches]
            if missing_jer:
                raise RuntimeError("Missing requested JER branch(es): " + ", ".join(missing_jer))
            branches.extend(jer_branches)
        missing_selection = [name for name in JET_SELECTION_BRANCHES if name not in all_branches]
        if missing_selection:
            raise RuntimeError(
                "DNN input construction requires selection branch(es): " + ", ".join(missing_selection)
            )
        branches.extend(JET_SELECTION_BRANCHES)
        weight_branches = [name for name in OPTIONAL_WEIGHT_BRANCHES if name in all_branches]
        branches.extend(weight_branches)

        branches = list(dict.fromkeys(branches))
        arrays = tree.arrays(branches, library="np")
        for source_name, alias_name in JET_BRANCH_ALIASES.items():
            arrays[alias_name] = arrays[source_name]
        arrays["Jet_btagScore_compat"] = arrays[jet_score_branch]
        arrays["Muon_promptMVA_compat"] = arrays[muon_prompt_branch]
        arrays["_jet_score_branch"] = jet_score_branch
        arrays["_muon_prompt_branch"] = muon_prompt_branch
        arrays["_muon_pt_branch"] = pt_branch
        arrays["_weight_branches"] = weight_branches
        gen_weight_sum = np.float32(1.0)
        if "genWeight" in weight_branches:
            if "genWeightSum" not in root_file:
                raise RuntimeError(f"Missing genWeightSum histogram in file: {file_name}")
            gen_weight_hist, _ = root_file["genWeightSum"].to_numpy()
            if len(gen_weight_hist) == 0 or float(gen_weight_hist[0]) == 0.0:
                raise RuntimeError(f"Invalid genWeightSum histogram in file: {file_name}")
            gen_weight_sum = np.float32(gen_weight_hist[0])
        arrays["_genWeightSum"] = gen_weight_sum
        return arrays


def resolve_default_btag_threshold(jet_score_branch):
    if jet_score_branch == "Jet_btagUParTAK4B":
        return 0.1272
    if jet_score_branch == "Jet_btagPNetB":
        return 0.1917
    raise RuntimeError(f"No default btag threshold for score branch: {jet_score_branch}")


def apply_jet_pt_variation(
    arrays,
    jet_pt_shift_branch=None,
    jet_pt_shift_direction="up",
    jer_direction=None,
    muon_pt_branch=None,
    jet_pt_threshold=30.0,
    btag_threshold=None,
    scale_jet_mass_with_pt=False,
):
    if jet_pt_shift_branch is not None and jer_direction is not None:
        raise RuntimeError("Use either --jet-pt-shift-branch or --jer-direction, not both.")

    nominal_pt = arrays["Jet_pt_JEC"]
    if jet_pt_shift_branch is not None:
        sign = 1.0 if jet_pt_shift_direction == "up" else -1.0
        varied_pt = nominal_pt * (1.0 + sign * arrays[jet_pt_shift_branch])
        variation_label = f"{jet_pt_shift_branch}_{jet_pt_shift_direction}"
    else:
        if jer_direction is not None:
            varied_pt = nominal_pt / arrays["Jet_JER_corr"] * arrays[f"Jet_JER_corr_{jer_direction}"]
            variation_label = f"JER_corr_{jer_direction}"
        else:
            varied_pt = nominal_pt
            variation_label = "nominal"

    if btag_threshold is None:
        btag_threshold = resolve_default_btag_threshold(arrays["_jet_score_branch"])

    arrays["Jet_pt_JEC"] = varied_pt
    arrays["Jet_pt"] = varied_pt
    if scale_jet_mass_with_pt:
        safe_ratio = varied_pt / np.maximum(nominal_pt, 1.0e-6)
        arrays["Jet_mass_JEC"] = arrays["Jet_mass_JEC"] * safe_ratio
        arrays["Jet_mass"] = arrays["Jet_mass_JEC"]

    good_masks = []
    bjet_masks = []
    for event_index in range(len(varied_pt)):
        good = (
            (np.asarray(varied_pt[event_index], dtype=np.float32) > jet_pt_threshold)
            & (np.asarray(arrays["Jet_rawFactor"][event_index], dtype=np.float32) < 0.9)
            & np.asarray(arrays["JetIdTight"][event_index], dtype=bool)
            & (np.asarray(arrays["Jet_drFromMuon"][event_index], dtype=np.float32) > 0.4)
        )
        bjet = good & (np.asarray(arrays["Jet_btagScore_compat"][event_index], dtype=np.float32) > btag_threshold)
        good_masks.append(good)
        bjet_masks.append(bjet)
    arrays["GoodJetCond"] = np.asarray(good_masks, dtype=object)
    arrays["BJetCond"] = np.asarray(bjet_masks, dtype=object)
    arrays["_jet_pt_variation"] = variation_label
    arrays["_jet_pt_threshold"] = np.float32(jet_pt_threshold)
    arrays["_btag_threshold"] = np.float32(btag_threshold)
    return arrays


def load_xs_lumi(sample_name, era, input_path):
    source_root = Path(__file__).resolve().parents[2]
    with open(source_root / "json" / "XS" / "Run3.json", encoding="utf-8") as handle:
        xs_map = json.load(handle)
    with open(source_root / "json" / "Lumi" / "Run3.json", encoding="utf-8") as handle:
        lumi_map = json.load(handle)

    resolved_sample_name = resolve_sample_name(sample_name, input_path, xs_map)
    if sample_name not in xs_map:
        sample_name = resolved_sample_name
    if sample_name not in xs_map:
        raise RuntimeError(f"Sample {sample_name} not found in XS json.")
    if era not in lumi_map:
        raise RuntimeError(f"Era {era} not found in lumi json.")
    return sample_name, np.float32(xs_map[sample_name]), np.float32(lumi_map[era])


def compute_weight(arrays, event_index, xs, lumi):
    # XS is stored in pb while lumi is in /fb, so include 1e3 to obtain yields.
    weight = float(xs) * float(lumi) * 1000.0
    for branch in arrays["_weight_branches"]:
        branch_value = float(arrays[branch][event_index])
        if branch == "genWeight":
            branch_value /= float(arrays["_genWeightSum"])
        weight *= branch_value
    return np.float32(weight)


def extract_event_row(arrays, event_index, event_weight):
    good_mask = np.asarray(arrays["GoodJetCond"][event_index], dtype=bool)
    bjet_mask = np.asarray(arrays["BJetCond"][event_index], dtype=bool)
    jet_pt = np.asarray(arrays["Jet_pt"][event_index], dtype=np.float32)
    jet_eta = np.asarray(arrays["Jet_eta"][event_index], dtype=np.float32)
    jet_phi = np.asarray(arrays["Jet_phi"][event_index], dtype=np.float32)
    jet_mass = np.asarray(arrays["Jet_mass"][event_index], dtype=np.float32)
    jet_score = np.asarray(arrays["Jet_btagScore_compat"][event_index], dtype=np.float32)

    if good_mask.sum() < 4:
        return None

    good_indices = np.flatnonzero(good_mask)
    good_scores = jet_score[good_indices]
    ordered_indices = good_indices[np.argsort(good_scores)[::-1][:4]]

    lead_mu_idx = int(arrays["leadingMuonIdx"][event_index])
    sublead_mu_idx = int(arrays["subleadingMuonIdx"][event_index])
    muon_pt = np.asarray(arrays[arrays["_muon_pt_branch"]][event_index], dtype=np.float32)
    muon_eta = np.asarray(arrays["Muon_eta"][event_index], dtype=np.float32)
    muon_phi = np.asarray(arrays["Muon_phi"][event_index], dtype=np.float32)
    muon_mass = np.asarray(arrays["Muon_mass"][event_index], dtype=np.float32)
    muon_charge = np.asarray(arrays["Muon_charge"][event_index], dtype=np.float32)
    muon_mini_iso = np.asarray(arrays["Muon_miniPFRelIso_all"][event_index], dtype=np.float32)
    muon_jet_rel_iso = np.asarray(arrays["Muon_jetRelIso"][event_index], dtype=np.float32)
    muon_jet_df = np.asarray(arrays["Muon_jetDF"][event_index], dtype=np.float32)
    muon_prompt_mva = np.asarray(arrays["Muon_promptMVA_compat"][event_index], dtype=np.float32)

    if min(lead_mu_idx, sublead_mu_idx) < 0:
        return None
    if max(lead_mu_idx, sublead_mu_idx) >= len(muon_pt):
        return None

    selected_pt = jet_pt[ordered_indices]
    selected_eta = jet_eta[ordered_indices]
    selected_phi = jet_phi[ordered_indices]
    selected_mass = jet_mass[ordered_indices]
    selected_score = jet_score[ordered_indices]
    selected_px, selected_py, selected_pz, _ = four_vector_from_pt_eta_phi_mass(
        selected_pt,
        selected_eta,
        selected_phi,
        selected_mass,
    )

    leading_mu = (
        muon_pt[lead_mu_idx],
        muon_eta[lead_mu_idx],
        muon_phi[lead_mu_idx],
        muon_mass[lead_mu_idx],
    )
    subleading_mu = (
        muon_pt[sublead_mu_idx],
        muon_eta[sublead_mu_idx],
        muon_phi[sublead_mu_idx],
        muon_mass[sublead_mu_idx],
    )

    lead_mu_dr = delta_r(selected_eta, selected_phi, leading_mu[1], leading_mu[2]).astype(np.float32)
    sublead_mu_dr = delta_r(selected_eta, selected_phi, subleading_mu[1], subleading_mu[2]).astype(np.float32)
    bjet_node = np.stack(
        [selected_px, selected_py, selected_pz, selected_mass, selected_score, lead_mu_dr, sublead_mu_dr],
        axis=1,
    ).reshape(-1)

    bjet_dr = np.array(
        [delta_r(selected_eta[i], selected_phi[i], selected_eta[j], selected_phi[j]) for i, j in PAIR_INDEX],
        dtype=np.float32,
    )
    bjet_invmass = np.array(
        [
            invariant_mass(
                selected_pt[i],
                selected_eta[i],
                selected_phi[i],
                selected_mass[i],
                selected_pt[j],
                selected_eta[j],
                selected_phi[j],
                selected_mass[j],
            )
            for i, j in PAIR_INDEX
        ],
        dtype=np.float32,
    )

    lead_px, lead_py, lead_pz, _ = four_vector_from_pt_eta_phi_mass(*leading_mu)
    sub_px, sub_py, sub_pz, _ = four_vector_from_pt_eta_phi_mass(*subleading_mu)
    muon_p3 = np.array([lead_px, lead_py, lead_pz, sub_px, sub_py, sub_pz], dtype=np.float32)
    muon_extra_features = np.array(
        [
            muon_mini_iso[lead_mu_idx],
            muon_mini_iso[sublead_mu_idx],
            muon_jet_rel_iso[lead_mu_idx],
            muon_jet_rel_iso[sublead_mu_idx],
            muon_jet_df[lead_mu_idx],
            muon_jet_df[sublead_mu_idx],
            muon_prompt_mva[lead_mu_idx],
            muon_prompt_mva[sublead_mu_idx],
        ],
        dtype=np.float32,
    )
    dimuon_dr = np.array(
        [delta_r(leading_mu[1], leading_mu[2], subleading_mu[1], subleading_mu[2])],
        dtype=np.float32,
    )
    dimuon_mass = np.array(
        [invariant_mass(*leading_mu, *subleading_mu)],
        dtype=np.float32,
    )

    n_good_jet = np.array([np.sum(good_mask)], dtype=np.float32)
    n_bjet = np.array([np.sum(good_mask & bjet_mask)], dtype=np.float32)
    sumjet_pt = np.array([np.sum(jet_pt[good_mask])], dtype=np.float32)
    sumbjet_pt = np.array([np.sum(jet_pt[good_mask & bjet_mask])], dtype=np.float32)
    sumjet_energy = np.sum(energy_from_pt_eta_mass(jet_pt[good_mask], jet_eta[good_mask], jet_mass[good_mask]))
    jet_cent = np.array([sumjet_pt[0] / sumjet_energy if sumjet_energy > 0 else 0.0], dtype=np.float32)
    weight = np.array([event_weight], dtype=np.float32)
    is_same_sign = np.array(
        [1.0 if (muon_charge[lead_mu_idx] * muon_charge[sublead_mu_idx]) > 0 else 0.0],
        dtype=np.float32,
    )

    return np.concatenate(
        [
            bjet_node,
            bjet_dr,
            bjet_invmass,
            muon_p3,
            muon_extra_features,
            dimuon_dr,
            dimuon_mass,
            n_good_jet,
            n_bjet,
            sumjet_pt,
            sumbjet_pt,
            jet_cent,
            weight,
            is_same_sign,
        ]
    ).astype(np.float32)


def to_structured_array(rows):
    dtype = [(field_name, np.float32) for field_name in ALL_FIELD_NAMES]
    structured = np.empty(len(rows), dtype=dtype)
    if len(rows) == 0:
        return structured
    stacked = np.stack(rows).astype(np.float32)
    for field_index, field_name in enumerate(ALL_FIELD_NAMES):
        structured[field_name] = stacked[:, field_index]
    return structured


def extract_file(
    input_path,
    output_path,
    sample_name=None,
    era=None,
    tree_name="Events",
    verbose=False,
    jet_pt_shift_branch=None,
    jet_pt_shift_direction="up",
    jer_direction=None,
    muon_pt_branch=None,
    jet_pt_threshold=30.0,
    btag_threshold=None,
    scale_jet_mass_with_pt=False,
):
    arrays = load_arrays(
        input_path,
        tree_name=tree_name,
        jet_pt_shift_branch=jet_pt_shift_branch,
        jer_direction=jer_direction,
        muon_pt_branch=muon_pt_branch,
    )
    arrays = apply_jet_pt_variation(
        arrays,
        jet_pt_shift_branch=jet_pt_shift_branch,
        jet_pt_shift_direction=jet_pt_shift_direction,
        jer_direction=jer_direction,
        jet_pt_threshold=jet_pt_threshold,
        btag_threshold=btag_threshold,
        scale_jet_mass_with_pt=scale_jet_mass_with_pt,
    )
    resolved_sample_name, xs, lumi = load_xs_lumi(sample_name, era, input_path)
    if verbose:
        print("signal extraction input:", input_path)
        print("resolved sample:", resolved_sample_name)
        print("XS:", float(xs))
        print("lumi:", float(lumi))
        print("genWeightSum:", float(arrays["_genWeightSum"]))
        print("jet score branch:", arrays["_jet_score_branch"])
        print("muon prompt branch:", arrays["_muon_prompt_branch"])
        print("muon pt branch:", arrays["_muon_pt_branch"])
        print("weight branches:", arrays["_weight_branches"])
        if "_jet_pt_variation" in arrays:
            print("jet pt variation:", arrays["_jet_pt_variation"])
            print("jet pt threshold:", float(arrays["_jet_pt_threshold"]))
            print("btag threshold:", float(arrays["_btag_threshold"]))

    rows = []
    n_input = len(arrays["Jet_pt"])
    for event_index in range(n_input):
        event_weight = compute_weight(arrays, event_index, xs, lumi)
        row = extract_event_row(arrays, event_index, event_weight)
        if row is not None:
            rows.append(row)

    if rows:
        output = to_structured_array(rows)
    else:
        output = to_structured_array([])
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.save(output_path, output)

    print(f"input events: {n_input}")
    print(f"saved events: {output.shape[0]}")
    print(f"feature dimension: {len(FEATURE_FIELD_NAMES)}")
    print(f"resolved sample: {resolved_sample_name}")
    print(f"saved output: {output_path}")


def main():
    args = parse_args()
    extract_file(
        input_path=args.input,
        output_path=args.output,
        sample_name=args.sample_name,
        era=args.era,
        tree_name=args.tree_name,
        verbose=bool(args.verbose),
        jet_pt_shift_branch=args.jet_pt_shift_branch,
        jet_pt_shift_direction=args.jet_pt_shift_direction,
        jer_direction=args.jer_direction,
        muon_pt_branch=args.muon_pt_branch,
        jet_pt_threshold=args.jet_pt_threshold,
        btag_threshold=args.btag_threshold,
        scale_jet_mass_with_pt=args.scale_jet_mass_with_pt,
    )


if __name__ == "__main__":
    main()
