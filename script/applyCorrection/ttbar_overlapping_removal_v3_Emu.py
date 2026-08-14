import ROOT
import os
import correctionlib


def processing(rdf, recordedModules, branchArray, era, ds=""):
    if "ttbar_overlapping_removal_v3_Emu.C" not in recordedModules:
        correctionlib.register_pyroot_binding()
        this_dir = os.path.dirname(os.path.abspath(__file__))
        ROOT.gInterpreter.AddIncludePath(this_dir)
        ROOT.gInterpreter.ProcessLine('#include "ttbar_overlapping_removal_v3_Emu.C"')
        recordedModules.append("ttbar_overlapping_removal_v3_Emu.C")

    ds_upper = ds.upper()
    is_ttbb = "TTBB" in ds_upper
    is_ttbar = ("TTBAR" in ds_upper) and not is_ttbb
    if not is_ttbb and not is_ttbar:
        return rdf, recordedModules, branchArray

    if "GenPart_isMuonForOverlap" not in branchArray:
        rdf = rdf.Define(
            "GenPart_isMuonForOverlap",
            "BuildLeptonForOverlapMaskV3Emu(GenPart_pdgId, GenPart_status, GenPart_statusFlags, GenPart_pt, GenPart_eta)",
        )
        branchArray.append("GenPart_isMuonForOverlap")

    if "GenPart_bHadOrigin" not in branchArray:
        rdf = rdf.Define(
            "GenPart_bHadOrigin",
            "GetBHadSourceV3Emu(GenPart_pdgId, GenPart_genPartIdxMother)",
        )
        branchArray.append("GenPart_bHadOrigin")

    if "GenJet_minMuonDr" not in branchArray:
        rdf = rdf.Define(
            "GenJet_minMuonDr",
            "GetGenJetMinLeptonDrV3Emu(GenJet_eta, GenJet_phi, GenPart_isMuonForOverlap, GenPart_eta, GenPart_phi)",
        )
        branchArray.append("GenJet_minMuonDr")

    if "GenJet_bHadIdx" not in branchArray:
        rdf = rdf.Define(
            "GenJet_bHadIdx",
            "GetGenJetBHadV3Emu(GenJet_hadronFlavour, GenJet_pt, GenJet_eta, GenJet_phi, GenPart_bHadOrigin, GenPart_eta, GenPart_phi, GenJet_minMuonDr)",
        )
        branchArray.append("GenJet_bHadIdx")

    if "GenJet_bHadSource" not in branchArray:
        rdf = rdf.Define(
            "GenJet_bHadSource",
            "GetBGenJetsSourceV3Emu(GenJet_bHadIdx, GenPart_bHadOrigin)",
        )
        branchArray.append("GenJet_bHadSource")

    if "GenJet_isFiducialBJet" not in branchArray:
        rdf = rdf.Define(
            "GenJet_isFiducialBJet",
            "(GenJet_hadronFlavour == 5) && (GenJet_pt > 20.0f) && (abs(GenJet_eta) < 2.5f)",
        )
        branchArray.append("GenJet_isFiducialBJet")

    if "nAdditionalFiducialBJet" not in branchArray:
        rdf = rdf.Define(
            "nAdditionalFiducialBJet",
            "Nonzero(GenJet_isFiducialBJet && (GenJet_bHadSource != 6) && (GenJet_bHadSource > 0)).size()",
        )
        branchArray.append("nAdditionalFiducialBJet")

    if "nTopFiducialBJet" not in branchArray:
        rdf = rdf.Define(
            "nTopFiducialBJet",
            "Nonzero(GenJet_isFiducialBJet && (GenJet_bHadSource == 6)).size()",
        )
        branchArray.append("nTopFiducialBJet")

    if "nMuonMaskedFiducialBJet" not in branchArray:
        rdf = rdf.Define(
            "nMuonMaskedFiducialBJet",
            "Nonzero(GenJet_isFiducialBJet && (GenJet_bHadSource == -13)).size()",
        )
        branchArray.append("nMuonMaskedFiducialBJet")

    if is_ttbb:
        rdf = rdf.Filter("(nAdditionalFiducialBJet >= 1)")
    else:
        rdf = rdf.Filter("(nAdditionalFiducialBJet == 0)")

    return rdf, recordedModules, branchArray
