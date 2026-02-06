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
            rdf = rdf.Define("Jet_pt_JEScorr", "get_JES_corr_pt(Jet_pt, Jet_rawFactor, Jet_eta, Jet_phi)")
        else:
            rdf = rdf.Define("Jet_pt_JEScorr", "get_JES_corr_pt_v12(Jet_pt, Jet_rawFactor, Jet_eta, Jet_phi)")
        branchArray.append("Jet_pt_JEScorr")
    if "Jet_pt_JEC" not in branchArray:
        rdf = rdf.Define("Jet_pt_JEC", "get_JER_corr_pt(Jet_pt_JEScorr, Jet_eta, Jet_phi, Jet_genJetIdx, GenJet_pt, GenJet_eta, GenJet_phi, Rho_fixedGridRhoFastjetAll, event )")
        branchArray.append("Jet_pt_JEC")

    return rdf, recordedModules, branchArray