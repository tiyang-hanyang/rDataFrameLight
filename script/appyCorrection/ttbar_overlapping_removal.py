import ROOT
import os
import correctionlib

def processing(rdf, recordedModules, branchArray, era, ds=""):
    if "ttbar_overlapping_removal.C" not in recordedModules:
        correctionlib.register_pyroot_binding()
        this_dir = os.path.dirname(os.path.abspath(__file__))
        ROOT.gInterpreter.AddIncludePath(this_dir)
        ROOT.gInterpreter.ProcessLine('#include "ttbar_overlapping_removal.C"')
        recordedModules.append("ttbar_overlapping_removal.C")

    # tracing back the bhadrons's source
    if "GenPart_bHadOrigin" not in branchArray:
        rdf = rdf.Define('GenPart_bHadOrigin', "GetBHadSource(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_statusFlags)")
        branchArray.append("GenPart_bHadOrigin")
    if "GenJet_bHadIdx" not in branchArray:
        # a method to get rid of GenJet_hadronFlavour
        rdf = rdf.Define('GenJet_bHadIdx', "GetGenJetBHad(GenJet_hadronFlavour, GenJet_eta, GenJet_phi, GenPart_bHadOrigin, GenPart_eta, GenPart_phi)")
        # rdf = rdf.Define('GenJet_bHadIdx', "GetGenJetBHad_withoutGM(GenJet_eta, GenJet_phi, GenPart_bHadOrigin, GenPart_eta, GenPart_phi)")
        branchArray.append("GenJet_bHadIdx")
    # should not look at reco jet for the removal, instead, should look at genjet
    # if "Jet_bHadSource" not in branchArray:
    #     rdf = rdf.Define('Jet_bHadSource', "GetBJetsSource(Jet_genJetIdx, GenJet_bHadIdx, GenPart_bHadOrigin)")
    #     branchArray.append("Jet_bHadSource")
    if "GenJet_bHadSource" not in branchArray:
        rdf = rdf.Define('GenJet_bHadSource', "GetBGenJetsSource(GenJet_bHadIdx, GenPart_bHadOrigin)")
        branchArray.append("GenJet_bHadSource")
    if "nAdditionalBJet" not in branchArray:
        # -1 means not a b-jet, otherwise, 6 is from t/tbar (source label the abs source)
        rdf = rdf.Define("nAdditionalBJet", "Nonzero(GenJet_bHadSource != 6 && GenJet_bHadSource != -1).size()")
        branchArray.append("nAdditionalBJet")
    rdf = rdf.Filter("(nAdditionalBJet<2)")

    # the parton matching case

    # if "nAdditionalBQuark" not in branchArray:
    #     rdf = rdf.Define("nAdditionalBQuark", "countAddB(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_statusFlags)")
    #     branchArray.append("nAdditionalBQuark")

    return rdf, recordedModules, branchArray