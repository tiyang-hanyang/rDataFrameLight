import ROOT
import os
import correctionlib


def processing(rdf, recordedModules, branchArray, era="", ds=""):
    if "BTagCorr_shape.C" not in recordedModules:
        correctionlib.register_pyroot_binding()
        this_dir = os.path.dirname(os.path.abspath(__file__))
        ROOT.gInterpreter.AddIncludePath(this_dir)
        ROOT.gInterpreter.ProcessLine('#include "BTagCorr_shape.C"')
        recordedModules.append("BTagCorr_shape.C")

    ROOT.gInterpreter.ProcessLine('BTagShape_init("'+era+'")')

    shape_syst_sources = [
        "lf",
        "lfstats1",
        "lfstats2",
        "hf",
        "hfstats1",
        "hfstats2",
        "cferr1",
        "cferr2",
    ]
    if "Jet_PNet_SF_medium" not in branchArray:
        rdf = rdf.Define("Jet_PNet_SF_medium", "get_SF_shape(\"central\", Jet_hadronFlavour, Jet_eta, Jet_pt_JEC, Jet_btagPNetB, Jet_passJetIdTightLepVeto)")
        branchArray.append("Jet_PNet_SF_medium")
    for source in shape_syst_sources:
        for direction in ["up", "down"]:
            sf_syst_branch = f"Jet_PNet_SF_medium_{source}_{direction}"
            if sf_syst_branch not in branchArray:
                rdf = rdf.Define(
                    sf_syst_branch,
                    f"get_SF_shape_flavour_syst(\"{source}\", \"{direction}\", Jet_hadronFlavour, Jet_eta, Jet_pt_JEC, Jet_btagPNetB, Jet_passJetIdTightLepVeto)"
                )
                branchArray.append(sf_syst_branch)

    def define_weight(branch_name, sf_branch):
        return rdf.Define(branch_name, f"compute_total_weight_shape({sf_branch}, GoodJetCond)")

    if "btag_weight" not in branchArray:
        rdf = define_weight("btag_weight", "Jet_PNet_SF_medium")
        branchArray.append("btag_weight")
    for source in shape_syst_sources:
        for direction in ["up", "down"]:
            weight_branch = f"btag_weight_{source}_{direction}"
            sf_syst_branch = f"Jet_PNet_SF_medium_{source}_{direction}"
            if weight_branch not in branchArray:
                rdf = define_weight(weight_branch, sf_syst_branch)
                branchArray.append(weight_branch)

    return rdf, recordedModules, branchArray
