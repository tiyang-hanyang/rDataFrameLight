import ROOT
import os
import correctionlib

# according to the tutorial https://gitlab.cern.ch/cms-analysis/jme/jerc-application-tutorial/-/blob/master/ApplyOnNanoAOD/
# 2024 data only has the L2 & L3 JES correction
# 2024 MC has L2 correction and JER correction
def processing(rdf, recordedModules, branchArray, era, ds=""):
    if "JEC_MC.C" not in recordedModules:
        correctionlib.register_pyroot_binding()
        this_dir = os.path.dirname(os.path.abspath(__file__))
        ROOT.gInterpreter.AddIncludePath(this_dir)
        ROOT.gInterpreter.ProcessLine('#include "JEC_MC.C"')
        recordedModules.append("JEC_MC.C")
    ROOT.gInterpreter.ProcessLine('JEC_MC_init("'+era+'")')

    if "Jet_pt_JEScorr" not in branchArray:
        if "2024" in era or era == "Run3Summer23BPixNanoAODv12":
            rdf = rdf.Define("Jet_pt_JEScorr", "get_JES_corr_pt(Jet_pt, Jet_pt, Jet_rawFactor, Jet_eta, Jet_phi)")
            rdf = rdf.Define("Jet_mass_JEScorr", "get_JES_corr_pt(Jet_mass, Jet_pt, Jet_rawFactor, Jet_eta, Jet_phi)")
        else:
            rdf = rdf.Define("Jet_pt_JEScorr", "get_JES_corr_pt_v12(Jet_pt, Jet_pt, Jet_rawFactor, Jet_eta, Jet_phi)")
            rdf = rdf.Define("Jet_mass_JEScorr", "get_JES_corr_pt_v12(Jet_mass, Jet_pt, Jet_rawFactor, Jet_eta, Jet_phi)")
        branchArray.append("Jet_pt_JEScorr")
        branchArray.append("Jet_mass_JEScorr")
    if "Jet_JER_corr" not in branchArray:
        rdf = rdf.Define("Jet_JER_corr", 'get_JER_corr(Jet_pt_JEScorr, Jet_eta, Jet_phi, Jet_genJetIdx, GenJet_pt, GenJet_eta, GenJet_phi, Rho_fixedGridRhoFastjetAll, event, "nom")')
        rdf = rdf.Define("Jet_JER_corr_up", 'get_JER_corr(Jet_pt_JEScorr, Jet_eta, Jet_phi, Jet_genJetIdx, GenJet_pt, GenJet_eta, GenJet_phi, Rho_fixedGridRhoFastjetAll, event, "up")')
        rdf = rdf.Define("Jet_JER_corr_down", 'get_JER_corr(Jet_pt_JEScorr, Jet_eta, Jet_phi, Jet_genJetIdx, GenJet_pt, GenJet_eta, GenJet_phi, Rho_fixedGridRhoFastjetAll, event, "down")')
        branchArray.append("Jet_JER_corr")
        branchArray.append("Jet_JER_corr_up")
        branchArray.append("Jet_JER_corr_down")
    if "Jet_pt_JEC" not in branchArray:
        rdf = rdf.Define("Jet_pt_JEC", "Jet_pt_JEScorr * Jet_JER_corr")
        rdf = rdf.Define("Jet_mass_JEC", "Jet_mass_JEScorr * Jet_JER_corr")
        branchArray.append("Jet_pt_JEC")
        branchArray.append("Jet_mass_JEC")

    jes_year = None
    if era in ["RunIII2024Summer24NanoAODv15", "RunIII2024Summer24NanoAODv15_SSCR", "RunIII2024Summer24NanoAODv15_AR"]:
        jes_year = "2024"
    elif era == "Run3Summer23BPixNanoAODv12":
        jes_year = "2023BPix"
    elif era == "Run3Summer23NanoAODv12":
        jes_year = "2023"
    elif era == "Run3Summer22EENanoAODv12":
        jes_year = "2022EE"
    elif era == "Run3Summer22NanoAODv12":
        jes_year = "2022"

    if jes_year is not None:
        jes_sources = [
            "CMS_scale_j_FlavorQCD",
            "CMS_scale_j_RelativeBal",
            "CMS_scale_j_HF",
            "CMS_scale_j_BBEC1",
            "CMS_scale_j_EC2",
            "CMS_scale_j_Absolute",
            f"CMS_scale_j_Absolute_{jes_year}",
            f"CMS_scale_j_HF_{jes_year}",
            f"CMS_scale_j_EC2_{jes_year}",
            f"CMS_scale_j_RelativeSample_{jes_year}",
            f"CMS_scale_j_BBEC1_{jes_year}",
        ]
        for source in jes_sources:
            if source not in branchArray:
                rdf = rdf.Define(source, f'get_JES_uncertainty("{source}", Jet_eta, Jet_pt_JEC)')
                branchArray.append(source)

    return rdf, recordedModules, branchArray
