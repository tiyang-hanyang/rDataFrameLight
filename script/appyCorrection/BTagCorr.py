import ROOT
import os
import correctionlib

# from BTV https://btv-wiki.docs.cern.ch/PerformanceCalibration/fixedWPSFRecommendations/
# we need to have the weight defined more carefully
# weight = \Product_btagged {SF_i*eff_i / eff_i} 
#       * \Product_nontagged {( 1 - SF_i*eff_i) / (1 - eff_i)} 
def processing(rdf, recordedModules, branchArray, era="RunIII2024Summer24NanoAODv15", ds=""):
    if "BTagCorr.C" not in recordedModules:
        correctionlib.register_pyroot_binding()
        this_dir = os.path.dirname(os.path.abspath(__file__))
        ROOT.gInterpreter.AddIncludePath(this_dir)
        ROOT.gInterpreter.ProcessLine('#include "BTagCorr.C"')
        ROOT.gInterpreter.ProcessLine('Btag_init("'+era+'")')
        recordedModules.append("BTagCorr.C")

    # the scale would depend on the particular DS, so need additional input appended, for other corrections, ds is not used and the default empty input is fine
    ROOT.gInterpreter.ProcessLine('reload_eff("'+era+'", "' + ds +'");')

    if "Jet_UParTB_SF_medium" not in branchArray:
        if "RunIII2024Summer24NanoAODv15" in era:
            rdf = rdf.Define("Jet_UParTB_SF_medium", "get_SF_fixedWP(\"central\", Jet_hadronFlavour, Jet_eta, Jet_pt_JEC)")
        else:
            rdf = rdf.Define("Jet_UParTB_SF_medium", "get_SF_shape(\"central\", Jet_hadronFlavour, Jet_eta, Jet_pt_JEC, Jet_btagPNetB)")
        branchArray.append("Jet_UParTB_SF_medium")
    if "Jet_btagEff_medium" not in branchArray:
        rdf = rdf.Define("Jet_btagEff_medium", "get_eff(\"central\", Jet_hadronFlavour, Jet_eta, Jet_pt_JEC)")
        branchArray.append("Jet_btagEff_medium")
    if "btag_weight" not in branchArray:
        if "RunIII2024Summer24NanoAODv15" in era:
            rdf = rdf.Define("btag_weight", "compute_total_weight_2024(Jet_UParTB_SF_medium, Jet_btagEff_medium, Jet_btagUParTAK4B, GoodJetCond)")
        elif era=="Run3Summer23NanoAODv12":
            rdf = rdf.Define("btag_weight", "compute_total_weight_2023(Jet_UParTB_SF_medium, Jet_btagEff_medium, Jet_btagPNetB, GoodJetCond)")
        elif era=="Run3Summer23BPixNanoAODv12":
            rdf = rdf.Define("btag_weight", "compute_total_weight_2023BPix(Jet_UParTB_SF_medium, Jet_btagEff_medium, Jet_btagPNetB, GoodJetCond)")
        elif era=="Run3Summer22NanoAODv12":
            rdf = rdf.Define("btag_weight", "compute_total_weight_2022(Jet_UParTB_SF_medium, Jet_btagEff_medium, Jet_btagPNetB, GoodJetCond)")
        elif era=="Run3Summer22EENanoAODv12":
            rdf = rdf.Define("btag_weight", "compute_total_weight_2022EE(Jet_UParTB_SF_medium, Jet_btagEff_medium, Jet_btagPNetB, GoodJetCond)")
        else:
            print("campaign not recognized:", era)
            exit()
        branchArray.append("btag_weight")

    return rdf, recordedModules, branchArray