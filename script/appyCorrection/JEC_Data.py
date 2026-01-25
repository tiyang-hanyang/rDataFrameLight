import ROOT
import os
import correctionlib

# according to the tutorial https://gitlab.cern.ch/cms-analysis/jme/jerc-application-tutorial/-/blob/master/ApplyOnNanoAOD/
# 2024 data only has the L2 & L3 JES correction
# 2024 MC has L2 correction and JER correction
def processing(rdf, recordedModules, branchArray, era, ds=""):
    if "JEC_Data.C" not in recordedModules:
        correctionlib.register_pyroot_binding()
        this_dir = os.path.dirname(os.path.abspath(__file__))
        ROOT.gInterpreter.AddIncludePath(this_dir)
        ROOT.gInterpreter.ProcessLine('#include "JEC_Data.C"')
        recordedModules.append("JEC_Data.C")
        ROOT.gInterpreter.ProcessLine('JEC_Data_init("'+era+'")')

    if "Jet_pt_JEC" not in branchArray:
        if "2024" in era or era == "Run2023D":
            rdf = rdf.Define("Jet_pt_JEC", "get_JES_corr_pt(run, Jet_pt, Jet_rawFactor, Jet_eta, Jet_phi)")
        elif era == "Run2023C":
            rdf = rdf.Define("Jet_pt_JEC", "get_JES_corr_pt_2023C(run, Jet_pt, Jet_rawFactor, Jet_eta, Jet_phi)")
        else:
            rdf = rdf.Define("Jet_pt_JEC", "get_JES_corr_pt_v12(run, Jet_pt, Jet_rawFactor, Jet_eta, Jet_phi)")
        branchArray.append("Jet_pt_JEC")

    return rdf, recordedModules, branchArray