import awkward as ak
import numpy as np
import uproot


REQUIRED_BRANCHES = [
    "Jet_pt_JEC",
    "Jet_eta",
    "Jet_phi",
    "Jet_mass_JEC",
    "Jet_neHEF",
    "Jet_neEmEF",
    "Jet_nConstituents",
    "Jet_muEF",
    "Jet_chHEF",
    "Jet_chMultiplicity",
    "Jet_chEmEF",
    "Jet_btagUParTAK4B",
    "Jet_btagPNetB",
    "Jets_bsource",
    "Muon_pt",
    "Muon_eta",
    "Muon_phi",
    "Muon_mass",
    "Muon_dxy",
    "Muon_dz",
    "Muon_sip3d",
    "Muon_jetDF",
    "Muon_mediumId",
    "Muon_miniPFRelIso_all",
    "Muon_promptMVA",
    "Muon_charge",
]

INFERENCE_REQUIRED_BRANCHES = [
    "run",
    "luminosityBlock",
    "event",
    "Jet_pt_JEC",
    "Jet_eta",
    "Jet_phi",
    "Jet_mass_JEC",
    "Jet_btagUParTAK4B",
    "BJetCond",
    "Muon_pt",
    "Muon_eta",
    "Muon_phi",
    "Muon_mass",
    "leadingMuonIdx",
    "subleadingMuonIdx",
]

JET_BRANCH_ALIASES = {
    "Jet_pt_JEC": "Jet_pt",
    "Jet_mass_JEC": "Jet_mass",
}


def akTocartesian(j):
    px = j.pt * np.cos(j.phi)
    py = j.pt * np.sin(j.phi)
    pz = j.pt * np.sinh(j.eta)
    energy = np.sqrt(px**2 + py**2 + pz**2 + j.mass**2)
    return px, py, pz, energy


def akGetInv(jet_pairs):
    j1 = jet_pairs.j1
    j2 = jet_pairs.j2
    px1, py1, pz1, energy1 = akTocartesian(j1)
    px2, py2, pz2, energy2 = akTocartesian(j2)
    px = px1 + px2
    py = py1 + py2
    pz = pz1 + pz2
    energy = energy1 + energy2
    m2 = energy**2 - (px**2 + py**2 + pz**2)
    m2 = ak.where(m2 > 0, m2, 0)
    return np.sqrt(m2) / 125.0


def akGetDR(jet_pairs):
    j1 = jet_pairs.j1
    j2 = jet_pairs.j2
    deta = j1.eta - j2.eta
    dphi = j1.phi - j2.phi
    dphi_abs = np.abs(dphi)
    dphi_wrapped = ak.where(dphi_abs > np.pi, 2 * np.pi - dphi_abs, dphi_abs)
    return np.sqrt(deta**2 + dphi_wrapped**2)


def akGetJetMuonInv(jets, muon):
    px, py, pz, energy = akTocartesian(jets)
    mu_px, mu_py, mu_pz, mu_energy = akTocartesian(muon)
    total_px = px + mu_px
    total_py = py + mu_py
    total_pz = pz + mu_pz
    total_energy = energy + mu_energy
    m2 = total_energy**2 - (total_px**2 + total_py**2 + total_pz**2)
    m2 = ak.where(m2 > 0, m2, 0)
    return np.sqrt(m2) / 125.0


def akGetJetMuonDR(jets, muon):
    deta = jets.eta - muon.eta
    dphi = jets.phi - muon.phi
    dphi_abs = np.abs(dphi)
    dphi_wrapped = ak.where(dphi_abs > np.pi, 2 * np.pi - dphi_abs, dphi_abs)
    return np.sqrt(deta**2 + dphi_wrapped**2)


def akJetAwayFromMuons(jets, muons):
    mu1 = muons[:, 0]
    mu2 = muons[:, 1]
    dR1 = akGetJetMuonDR(jets, mu1)
    dR2 = akGetJetMuonDR(jets, mu2)
    return (dR1 > 0.4) & (dR2 > 0.4)


def akMuonSel(muons):
    return (
        (muons.pt > 15.0)
        & (abs(muons.eta) < 2.4)
        & (abs(muons.dxy) < 0.05)
        & (abs(muons.dz) < 0.1)
        & (abs(muons.sip3d) < 8)
        & (muons.mediumId == 1)
        & (muons.miniIso < 0.4)
        & (muons.jetDF < 0.2480)
        & (muons.promptMVA > 0.64)
    )


def akJetSel(jets):
    return (
        (jets.pt > 30.0)
        & (abs(jets.eta) < 2.5)
        & (jets.neHEF < 0.99)
        & (jets.neEmEF < 0.90)
        & (jets.nConstituents > 1)
        & (jets.muEF < 0.80)
        & (jets.chHEF > 0.01)
        & (jets.chMultiplicity > 0)
        & (jets.chEmEF < 0.80)
        & (jets.neHEF + jets.chEmEF < 0.90)
    )


def akbtagSel(jets):
    return jets.btagUParTAK4B > 0.1272


def require_branches(tree, required_branches):
    all_branches = {branch.name for branch in tree.branches}
    missing = [name for name in required_branches if name not in all_branches]
    if missing:
        raise RuntimeError("Missing required branches: " + ", ".join(missing))


def load_event_arrays(file_name, branches=None, tree_name="Events"):
    selected_branches = REQUIRED_BRANCHES if branches is None else branches
    with uproot.open(file_name) as root_file:
        tree = root_file[tree_name]
        require_branches(tree, selected_branches)
        arrays = tree.arrays(selected_branches, library="ak")
        for source_name, alias_name in JET_BRANCH_ALIASES.items():
            if source_name in arrays.fields:
                arrays = ak.with_field(arrays, arrays[source_name], alias_name)
        return arrays


def build_muons(info):
    return ak.zip(
        {
            "pt": info["Muon_pt"],
            "eta": info["Muon_eta"],
            "phi": info["Muon_phi"],
            "mass": info["Muon_mass"],
            "dxy": info["Muon_dxy"],
            "dz": info["Muon_dz"],
            "sip3d": info["Muon_sip3d"],
            "jetDF": info["Muon_jetDF"],
            "mediumId": info["Muon_mediumId"],
            "miniIso": info["Muon_miniPFRelIso_all"],
            "promptMVA": info["Muon_promptMVA"],
            "charge": info["Muon_charge"],
        },
        with_name="Muon",
    )


def build_jets(info):
    return ak.zip(
        {
            "pt": info["Jet_pt"],
            "eta": info["Jet_eta"],
            "phi": info["Jet_phi"],
            "mass": info["Jet_mass"],
            "neHEF": info["Jet_neHEF"],
            "neEmEF": info["Jet_neEmEF"],
            "nConstituents": info["Jet_nConstituents"],
            "muEF": info["Jet_muEF"],
            "chHEF": info["Jet_chHEF"],
            "chMultiplicity": info["Jet_chMultiplicity"],
            "chEmEF": info["Jet_chEmEF"],
            "btagUParTAK4B": info["Jet_btagUParTAK4B"],
            "source": info["Jets_bsource"],
            "PNetScore": info["Jet_btagPNetB"],
        },
        with_name="Jet",
    )


def build_inference_jets(info):
    return ak.zip(
        {
            "pt": info["Jet_pt"],
            "eta": info["Jet_eta"],
            "phi": info["Jet_phi"],
            "mass": info["Jet_mass"],
            "btagUParTAK4B": info["Jet_btagUParTAK4B"],
            "isSelected": info["BJetCond"],
        },
        with_name="Jet",
    )


def build_inference_muons(info):
    return ak.zip(
        {
            "pt": info["Muon_pt"],
            "eta": info["Muon_eta"],
            "phi": info["Muon_phi"],
            "mass": info["Muon_mass"],
        },
        with_name="Muon",
    )


def select_dimuon_events(info, muons):
    is_good_muon = akMuonSel(muons)
    n_good_muon = ak.sum(is_good_muon, axis=1)
    dimuon_mask = n_good_muon > 1

    dimuon_info = info[dimuon_mask]
    dimuon_muons = muons[dimuon_mask]
    dimuon_good_muon = is_good_muon[dimuon_mask]
    n_good_muon = n_good_muon[dimuon_mask]

    good_charge = dimuon_muons.charge[dimuon_good_muon]
    sum_charge = ak.sum(good_charge[:, :2], axis=1)
    pass_ss = (n_good_muon == 2) & (np.abs(sum_charge) == 2)

    return {
        "info": dimuon_info,
        "muons": dimuon_muons,
        "nGoodMuon": n_good_muon,
        "passSS": pass_ss,
    }


def select_bjet_events(selection):
    jets = build_jets(selection["info"])
    jets = jets[akJetSel(jets)]
    jets = jets[akJetAwayFromMuons(jets, selection["muons"])]

    n_good_jet = ak.num(jets)
    event_mask = n_good_jet > 3

    jets = jets[event_mask]
    n_good_jet = n_good_jet[event_mask]
    muons = selection["muons"][event_mask]
    n_good_muon = selection["nGoodMuon"][event_mask]
    pass_ss = selection["passSS"][event_mask]

    total_good_jet_pt = ak.sum(jets.pt, axis=1)
    order = ak.argsort(jets.btagUParTAK4B, axis=1, ascending=False)
    jets = jets[order]

    is_bjet = akbtagSel(jets)
    n_bjet = ak.sum(is_bjet, axis=1)
    total_bjet_pt = ak.sum(jets.pt[is_bjet], axis=1)

    return {
        "jets": jets[:, :4],
        "muons": muons,
        "nGoodJet": n_good_jet,
        "nBJet": n_bjet,
        "totalGoodJetPt": total_good_jet_pt,
        "totalBJetPt": total_bjet_pt,
        "nGoodMuon": n_good_muon,
        "passSS": pass_ss,
    }


def compute_bjet_assignment_labels(jet_pairs):
    both_higgs = (jet_pairs.j1.source == 25) & (jet_pairs.j2.source == 25)
    matched = ak.any(both_higgs, axis=1)
    first_true = ak.argmax(both_higgs, axis=1)
    return ak.where(matched, first_true, -1)


def compute_feature_blocks(selected_events):
    jets = selected_events["jets"]
    muons = selected_events["muons"]
    jet_pairs = ak.combinations(jets, 2, axis=1, fields=["j1", "j2"])

    inv_mass = akGetInv(jet_pairs)
    delta_r = akGetDR(jet_pairs)

    leading_muon = muons[:, 0]
    subleading_muon = muons[:, 1]

    jet_leading_inv = akGetJetMuonInv(jets, leading_muon)
    jet_leading_dr = akGetJetMuonDR(jets, leading_muon)
    jet_subleading_inv = akGetJetMuonInv(jets, subleading_muon)
    jet_subleading_dr = akGetJetMuonDR(jets, subleading_muon)

    return {
        "jets": jets,
        "invMass": inv_mass,
        "deltaR": delta_r,
        "jetLeadingInv": jet_leading_inv,
        "jetLeadingDR": jet_leading_dr,
        "jetSubleadingInv": jet_subleading_inv,
        "jetSubleadingDR": jet_subleading_dr,
        "bjetCat": compute_bjet_assignment_labels(jet_pairs),
        "nGoodJet": selected_events["nGoodJet"],
        "nBJet": selected_events["nBJet"],
        "totalGoodJetPt": selected_events["totalGoodJetPt"],
        "totalBJetPt": selected_events["totalBJetPt"],
        "nGoodMuon": selected_events["nGoodMuon"],
        "passSS": selected_events["passSS"],
    }


def select_inference_events(info):
    jets = build_inference_jets(info)
    muons = build_inference_muons(info)

    selected_jets = jets[jets.isSelected]
    order = ak.argsort(selected_jets.btagUParTAK4B, axis=1, ascending=False)
    selected_jets = selected_jets[order]

    event_mask = ak.num(selected_jets) > 3
    selected_jets = selected_jets[event_mask][:, :4]
    selected_muons = muons[event_mask]
    selected_info = info[event_mask]

    leading_idx = selected_info["leadingMuonIdx"]
    subleading_idx = selected_info["subleadingMuonIdx"]
    event_index = ak.local_index(leading_idx, axis=0)
    leading_muons = selected_muons[event_index, leading_idx]
    subleading_muons = selected_muons[event_index, subleading_idx]

    return {
        "ids": {
            "run": selected_info["run"],
            "luminosityBlock": selected_info["luminosityBlock"],
            "event": selected_info["event"],
        },
        "jets": selected_jets,
        "leadingMuon": leading_muons,
        "subleadingMuon": subleading_muons,
        "nGoodJet": ak.num(selected_jets),
        "nBJet": ak.sum(selected_jets.isSelected, axis=1),
        "totalGoodJetPt": ak.sum(selected_jets.pt, axis=1),
        "totalBJetPt": ak.sum(selected_jets.pt, axis=1),
    }


def compute_inference_feature_blocks(selected_events):
    jets = selected_events["jets"]
    jet_pairs = ak.combinations(jets, 2, axis=1, fields=["j1", "j2"])

    inv_mass = akGetInv(jet_pairs)
    delta_r = akGetDR(jet_pairs)
    jet_leading_inv = akGetJetMuonInv(jets, selected_events["leadingMuon"])
    jet_leading_dr = akGetJetMuonDR(jets, selected_events["leadingMuon"])
    jet_subleading_inv = akGetJetMuonInv(jets, selected_events["subleadingMuon"])
    jet_subleading_dr = akGetJetMuonDR(jets, selected_events["subleadingMuon"])

    return {
        "ids": selected_events["ids"],
        "jets": jets,
        "invMass": inv_mass,
        "deltaR": delta_r,
        "jetLeadingInv": jet_leading_inv,
        "jetLeadingDR": jet_leading_dr,
        "jetSubleadingInv": jet_subleading_inv,
        "jetSubleadingDR": jet_subleading_dr,
        "nGoodJet": selected_events["nGoodJet"],
        "nBJet": selected_events["nBJet"],
        "totalGoodJetPt": selected_events["totalGoodJetPt"],
        "totalBJetPt": selected_events["totalBJetPt"],
    }


def to_numpy_column(array):
    return ak.to_numpy(array)


def build_training_matrix(feature_blocks):
    jets = feature_blocks["jets"]
    jet_leading_inv = feature_blocks["jetLeadingInv"]
    jet_leading_dr = feature_blocks["jetLeadingDR"]
    jet_subleading_inv = feature_blocks["jetSubleadingInv"]
    jet_subleading_dr = feature_blocks["jetSubleadingDR"]
    delta_r = feature_blocks["deltaR"]
    inv_mass = feature_blocks["invMass"]

    columns = []
    for jet_index in range(4):
        columns.extend(
            [
                to_numpy_column(jets[:, jet_index].pt),
                to_numpy_column(jets[:, jet_index].eta),
                to_numpy_column(jets[:, jet_index].btagUParTAK4B),
                to_numpy_column(jet_leading_inv[:, jet_index]),
                to_numpy_column(jet_leading_dr[:, jet_index]),
                to_numpy_column(jet_subleading_inv[:, jet_index]),
                to_numpy_column(jet_subleading_dr[:, jet_index]),
            ]
        )

    columns.extend(
        [
            to_numpy_column(feature_blocks["nGoodJet"]),
            to_numpy_column(feature_blocks["nBJet"]),
            to_numpy_column(feature_blocks["totalGoodJetPt"]),
            to_numpy_column(feature_blocks["totalBJetPt"]),
        ]
    )

    for pair_index in range(6):
        columns.append(to_numpy_column(delta_r[:, pair_index]))
    for pair_index in range(6):
        columns.append(to_numpy_column(inv_mass[:, pair_index]))

    columns.extend(
        [
            to_numpy_column(feature_blocks["bjetCat"]),
            to_numpy_column(feature_blocks["nGoodMuon"]),
            to_numpy_column(feature_blocks["passSS"]),
        ]
    )

    return np.column_stack(columns)


def build_inference_matrix(feature_blocks, include_ids=False):
    jets = feature_blocks["jets"]
    jet_leading_inv = feature_blocks["jetLeadingInv"]
    jet_leading_dr = feature_blocks["jetLeadingDR"]
    jet_subleading_inv = feature_blocks["jetSubleadingInv"]
    jet_subleading_dr = feature_blocks["jetSubleadingDR"]
    delta_r = feature_blocks["deltaR"]
    inv_mass = feature_blocks["invMass"]

    columns = []
    for jet_index in range(4):
        columns.extend(
            [
                to_numpy_column(jets[:, jet_index].pt),
                to_numpy_column(jets[:, jet_index].eta),
                to_numpy_column(jets[:, jet_index].btagUParTAK4B),
                to_numpy_column(jet_leading_inv[:, jet_index]),
                to_numpy_column(jet_leading_dr[:, jet_index]),
                to_numpy_column(jet_subleading_inv[:, jet_index]),
                to_numpy_column(jet_subleading_dr[:, jet_index]),
            ]
        )

    columns.extend(
        [
            to_numpy_column(feature_blocks["nGoodJet"]),
            to_numpy_column(feature_blocks["nBJet"]),
            to_numpy_column(feature_blocks["totalGoodJetPt"]),
            to_numpy_column(feature_blocks["totalBJetPt"]),
        ]
    )

    for pair_index in range(6):
        columns.append(to_numpy_column(delta_r[:, pair_index]))
    for pair_index in range(6):
        columns.append(to_numpy_column(inv_mass[:, pair_index]))

    features = np.column_stack(columns)
    if not include_ids:
        return features

    return {
        "run": to_numpy_column(feature_blocks["ids"]["run"]),
        "luminosityBlock": to_numpy_column(feature_blocks["ids"]["luminosityBlock"]),
        "event": to_numpy_column(feature_blocks["ids"]["event"]),
        "features": features,
    }


def extract_features_from_arrays(info):
    muons = build_muons(info)
    dimuon_events = select_dimuon_events(info, muons)
    selected_events = select_bjet_events(dimuon_events)
    feature_blocks = compute_feature_blocks(selected_events)
    return build_training_matrix(feature_blocks)


def extract_features_from_file(file_name, tree_name="Events"):
    info = load_event_arrays(file_name, tree_name=tree_name)
    return extract_features_from_arrays(info)


def extract_inference_feature_blocks_from_arrays(info):
    selected_events = select_inference_events(info)
    return compute_inference_feature_blocks(selected_events)


def extract_inference_feature_blocks_from_file(file_name, tree_name="Events"):
    info = load_event_arrays(file_name, branches=INFERENCE_REQUIRED_BRANCHES, tree_name=tree_name)
    return extract_inference_feature_blocks_from_arrays(info)


def extract_inference_features_from_arrays(info, include_ids=False):
    feature_blocks = extract_inference_feature_blocks_from_arrays(info)
    return build_inference_matrix(feature_blocks, include_ids=include_ids)


def extract_inference_features_from_file(file_name, tree_name="Events", include_ids=False):
    feature_blocks = extract_inference_feature_blocks_from_file(file_name, tree_name=tree_name)
    return build_inference_matrix(feature_blocks, include_ids=include_ids)
