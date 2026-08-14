import ROOT
import os
import correctionlib


def processing(rdf, recordedModules, branchArray, era="RunIII2024Summer24NanoAODv15", ds=""):
    if "BTagCorr_fixedWP.C" not in recordedModules:
        correctionlib.register_pyroot_binding()
        this_dir = os.path.dirname(os.path.abspath(__file__))
        ROOT.gInterpreter.AddIncludePath(this_dir)
        ROOT.gInterpreter.ProcessLine('#include "BTagCorr_fixedWP.C"')
        recordedModules.append("BTagCorr_fixedWP.C")

    ROOT.gInterpreter.ProcessLine('BtagFixedWP_init("'+era+'", "' + ds +'");')

    wp_map = {
        "loose": "L",
        "medium": "M",
        "tight": "T",
    }
    wp_syst_sources = [
        "mur",
        "pileup",
        "pdfas",
        "isrdef",
        "topmass",
        "bfragmentation",
        "hdamp",
        "jer",
        "statistic",
        "muf",
        "fsrdef",
        "jes",
        "type3",
    ]
    for label, wp in wp_map.items():
        sf_branch = f"Jet_UParTB_SF_{label}"
        sf_branch_up = f"{sf_branch}_up"
        sf_branch_down = f"{sf_branch}_down"
        eff_branch = f"Jet_btagEff_{label}"
        weight_branch = f"btag_weight_{label}"
        weight_branch_up = f"{weight_branch}_up"
        weight_branch_down = f"{weight_branch}_down"
        threshold = {"loose": "0.0246", "medium": "0.1272", "tight": "0.4648"}[label]

        if sf_branch not in branchArray:
            rdf = rdf.Define(sf_branch, f'get_SF_fixedWP("central", "{wp}", Jet_hadronFlavour, Jet_eta, Jet_pt_JEC)')
            branchArray.append(sf_branch)
        if sf_branch_up not in branchArray:
            rdf = rdf.Define(sf_branch_up, f'get_SF_fixedWP("up", "{wp}", Jet_hadronFlavour, Jet_eta, Jet_pt_JEC)')
            branchArray.append(sf_branch_up)
        if sf_branch_down not in branchArray:
            rdf = rdf.Define(sf_branch_down, f'get_SF_fixedWP("down", "{wp}", Jet_hadronFlavour, Jet_eta, Jet_pt_JEC)')
            branchArray.append(sf_branch_down)
        if eff_branch not in branchArray:
            rdf = rdf.Define(eff_branch, f'get_eff("central", "{wp}", Jet_hadronFlavour, Jet_eta, Jet_pt_JEC)')
            branchArray.append(eff_branch)
        if weight_branch not in branchArray:
            rdf = rdf.Define(weight_branch, f"compute_total_weight_fixedWP_2024({sf_branch}, {eff_branch}, Jet_btagUParTAK4B, GoodJetCond, {threshold})")
            branchArray.append(weight_branch)
        if weight_branch_up not in branchArray:
            rdf = rdf.Define(weight_branch_up, f"compute_total_weight_fixedWP_2024({sf_branch_up}, {eff_branch}, Jet_btagUParTAK4B, GoodJetCond, {threshold})")
            branchArray.append(weight_branch_up)
        if weight_branch_down not in branchArray:
            rdf = rdf.Define(weight_branch_down, f"compute_total_weight_fixedWP_2024({sf_branch_down}, {eff_branch}, Jet_btagUParTAK4B, GoodJetCond, {threshold})")
            branchArray.append(weight_branch_down)

        if label == "medium":
            for source in wp_syst_sources:
                for direction in ["up", "down"]:
                    weight_syst_branch = f"btag_weight_{source}_{direction}"
                    if weight_syst_branch not in branchArray:
                        rdf = rdf.Define(
                            weight_syst_branch,
                            f'compute_total_weight_fixedWP_2024_syst("{direction}_{source}", "{wp}", Jet_hadronFlavour, Jet_eta, Jet_pt_JEC, Jet_btagUParTAK4B, GoodJetCond, {threshold})'
                        )
                        branchArray.append(weight_syst_branch)

    if "btag_weight" not in branchArray:
        rdf = rdf.Define("btag_weight", "btag_weight_medium")
        branchArray.append("btag_weight")
    if "btag_weight_up" not in branchArray:
        syst_weight_list = ", ".join([f"btag_weight_{source}_{direction}" for source in wp_syst_sources for direction in ["up", "down"]])
        rdf = rdf.Define("btag_weight_up", f"btagWeightQuadratureUp(btag_weight, std::vector<float>{{{syst_weight_list}}})")
        branchArray.append("btag_weight_up")
    if "btag_weight_down" not in branchArray:
        syst_weight_list = ", ".join([f"btag_weight_{source}_{direction}" for source in wp_syst_sources for direction in ["up", "down"]])
        rdf = rdf.Define("btag_weight_down", f"btagWeightQuadratureDown(btag_weight, std::vector<float>{{{syst_weight_list}}})")
        branchArray.append("btag_weight_down")

    return rdf, recordedModules, branchArray
