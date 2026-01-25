import ROOT
import os
import correctionlib

# muon scale correction only now, later consider about the corrections for the trigger
def processing(rdf, recordedModules, branchArray, era, ds=""):
    if "Muon_corr.C" not in recordedModules:
        correctionlib.register_pyroot_binding()
        this_dir = os.path.dirname(os.path.abspath(__file__))
        ROOT.gInterpreter.AddIncludePath(this_dir)
        ROOT.gInterpreter.ProcessLine('#include "Muon_corr.C"')
        recordedModules.append("Muon_corr.C")
        ROOT.gInterpreter.ProcessLine('Muon_corr_init("'+era+'")')

    if "Muon_IDscale" not in branchArray:
        rdf = rdf.Define("Muon_IDscale", "MuonIDScale(Muon_eta, Muon_pt_Rcorr)")
        branchArray.append("Muon_IDscale")
    if "Muon_Isoscale" not in branchArray:
        rdf = rdf.Define("Muon_Isoscale", "MuonIsoScale(Muon_eta, Muon_pt_Rcorr)")
        branchArray.append("Muon_Isoscale")
    if "Muon_MVAscale" not in branchArray:
        rdf = rdf.Define("Muon_MVAscale", "MuonMVAScale(Muon_eta, Muon_pt_Rcorr)")
        branchArray.append("Muon_MVAscale")
    if "MuonScale" not in branchArray:
        rdf = rdf.Define("MuonScale", "EventMuonSF(Muon_IDscale, Muon_Isoscale, Muon_MVAscale, leadingMuonIdx, subleadingMuonIdx)")
        branchArray.append("MuonScale")

    return rdf, recordedModules, branchArray