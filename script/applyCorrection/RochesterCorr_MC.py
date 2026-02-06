import ROOT
import os
import correctionlib

def processing(rdf, recordedModules, branchArray, era, ds=""):
    if "RochesterCorr_MC.C" not in recordedModules:
        correctionlib.register_pyroot_binding()
        this_dir = os.path.dirname(os.path.abspath(__file__))
        ROOT.gInterpreter.AddIncludePath(this_dir)
        ROOT.gInterpreter.ProcessLine('#include "RochesterCorr_MC.C"')
        recordedModules.append("RochesterCorr_MC.C")
        ROOT.gInterpreter.ProcessLine('RCorr_MC_init("'+era+'")')


    # scale correction
    if "Muon_pt_Rscale" not in branchArray:
        rdf = rdf.Define('Muon_pt_Rscale', 'pt_scale(0, Muon_pt, Muon_eta, Muon_phi, Muon_charge)')
        branchArray.append("Muon_pt_Rscale")
    # resolution correction
    if "Muon_pt_Rcorr" not in branchArray:
        rdf = rdf.Define('Muon_pt_Rcorr','pt_resol(Muon_pt_Rscale, Muon_eta, Muon_phi, Muon_nTrackerLayers, event, luminosityBlock)')
        branchArray.append("Muon_pt_Rcorr")
    # uncertainties
    if "Muon_pt_Rscale_up" not in branchArray:
        rdf = rdf.Define('Muon_pt_Rscale_up', 'pt_scale_var(Muon_pt_Rcorr, Muon_eta, Muon_phi, Muon_charge, "up")')
        branchArray.append("Muon_pt_Rscale_up")
    if "Muon_pt_Rscale_dn" not in branchArray:
        rdf = rdf.Define('Muon_pt_Rscale_dn', 'pt_scale_var(Muon_pt_Rcorr, Muon_eta, Muon_phi, Muon_charge, "dn")')
        branchArray.append("Muon_pt_Rscale_dn")
    if "Muon_pt_Rcorr_resolup" not in branchArray:
        rdf = rdf.Define("Muon_pt_Rcorr_resolup", 'pt_resol_var(Muon_pt_Rscale, Muon_pt_Rcorr, Muon_eta, "up")')
        branchArray.append("Muon_pt_Rcorr_resolup")
    if "Muon_pt_Rcorr_resoldn" not in branchArray:
        rdf = rdf.Define("Muon_pt_Rcorr_resoldn", 'pt_resol_var(Muon_pt_Rscale, Muon_pt_Rcorr, Muon_eta, "dn")')
        branchArray.append("Muon_pt_Rcorr_resoldn")

    return rdf, recordedModules, branchArray