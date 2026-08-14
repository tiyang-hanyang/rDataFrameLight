import ROOT

JES_SOURCES = [
    "CMS_scale_j_FlavorQCD",
    "CMS_scale_j_RelativeBal",
    "CMS_scale_j_HF",
    "CMS_scale_j_BBEC1",
    "CMS_scale_j_EC2",
    "CMS_scale_j_Absolute",
    "CMS_scale_j_Absolute_2022",
    "CMS_scale_j_HF_2022",
    "CMS_scale_j_EC2_2022",
    "CMS_scale_j_RelativeSample_2022",
    "CMS_scale_j_BBEC1_2022",
    "CMS_scale_j_Absolute_2022EE",
    "CMS_scale_j_HF_2022EE",
    "CMS_scale_j_EC2_2022EE",
    "CMS_scale_j_RelativeSample_2022EE",
    "CMS_scale_j_BBEC1_2022EE",
    "CMS_scale_j_Absolute_2023",
    "CMS_scale_j_HF_2023",
    "CMS_scale_j_EC2_2023",
    "CMS_scale_j_RelativeSample_2023",
    "CMS_scale_j_BBEC1_2023",
    "CMS_scale_j_Absolute_2023BPix",
    "CMS_scale_j_HF_2023BPix",
    "CMS_scale_j_EC2_2023BPix",
    "CMS_scale_j_RelativeSample_2023BPix",
    "CMS_scale_j_BBEC1_2023BPix",
    "CMS_scale_j_Absolute_2024",
    "CMS_scale_j_HF_2024",
    "CMS_scale_j_EC2_2024",
    "CMS_scale_j_RelativeSample_2024",
    "CMS_scale_j_BBEC1_2024",
]


def _append_once(branchArray, branch):
    if branch not in branchArray:
        branchArray.append(branch)


def _uses_v15_jetid(era):
    # 2022/2023 samples now use the v15-style NanoAOD jet ID interface.
    return ("2022" in era) or ("2023" in era) or ("2024" in era) or ("NanoAODv15" in era)


def _define_jvm_weight_variation(rdf, branchArray, suffix, jet_pt_expr):
    loose_branch = f"looseJetCond_{suffix}"
    fail_branch = f"isLooseJetFailJVM_{suffix}"
    weight_branch = f"JVMweight_{suffix}"
    if weight_branch in branchArray:
        return rdf, branchArray
    rdf = rdf.Define(loose_branch, f"({jet_pt_expr} > 15.0) && JetIdTight && (Jet_neEmEF + Jet_chEmEF < 0.9)")
    rdf = rdf.Define(fail_branch, f"({loose_branch} && !(PassJetVeto))")
    rdf = rdf.Define(weight_branch, f"float(Nonzero({fail_branch}).size()==0)")
    _append_once(branchArray, weight_branch)
    return rdf, branchArray


def processing(rdf, recordedModules, branchArray, era, ds=""):
    if "JetVetoMap_eMu.C" not in recordedModules:
        ROOT.gInterpreter.ProcessLine('#include "JetVetoMap_eMu.C"')
        recordedModules.append("JetVetoMap_eMu.C")
    ROOT.gInterpreter.ProcessLine('JVM_eMu_init("' + era + '")')

    if "PassJetVeto" not in branchArray:
        rdf = rdf.Define("PassJetVeto", "passJetVetoFunc_eMu(Jet_eta, Jet_phi)")
        _append_once(branchArray, "PassJetVeto")
    if "Jet_passJetIdTight" not in branchArray:
        if _uses_v15_jetid(era):
            rdf = rdf.Define(
                "Jet_passJetIdTight",
                """
                ((abs(Jet_eta) <= 2.6) && (Jet_neHEF < 0.99) && (Jet_neEmEF < 0.90) && ((Jet_chMultiplicity + Jet_neMultiplicity) > 1) && (Jet_chHEF > 0.01) && (Jet_chMultiplicity > 0))
                || ((abs(Jet_eta) > 2.6) && (abs(Jet_eta) <= 2.7) && (Jet_neHEF < 0.90) && (Jet_neEmEF < 0.99))
                || ((abs(Jet_eta) > 2.7) && (abs(Jet_eta) <= 3.0) && (Jet_neHEF < 0.99))
                || ((abs(Jet_eta) > 3.0) && (Jet_neMultiplicity >= 2) && (Jet_neEmEF < 0.4))
                """,
            )
        else:
            rdf = rdf.Define(
                "Jet_passJetIdTight",
                """
                ((abs(Jet_eta) <= 2.7) && ((Jet_jetId & (1 << 1)) != 0))
                || ((abs(Jet_eta) > 2.7) && (abs(Jet_eta) <= 3.0) && ((Jet_jetId & (1 << 1)) != 0) && (Jet_neHEF < 0.99))
                || ((abs(Jet_eta) > 3.0) && ((Jet_jetId & (1 << 1)) != 0) && (Jet_neEmEF < 0.4))
                """,
            )
        _append_once(branchArray, "Jet_passJetIdTight")
    if "Jet_passJetIdTightLepVeto" not in branchArray:
        rdf = rdf.Define(
            "Jet_passJetIdTightLepVeto",
            """
            ((abs(Jet_eta) <= 2.7) && Jet_passJetIdTight && (Jet_muEF < 0.8) && (Jet_chEmEF < 0.8))
            || ((abs(Jet_eta) > 2.7) && Jet_passJetIdTight)
            """,
        )
        _append_once(branchArray, "Jet_passJetIdTightLepVeto")
    if "JetIdTight" not in branchArray:
        rdf = rdf.Define("JetIdTight", "Jet_passJetIdTightLepVeto")
        _append_once(branchArray, "JetIdTight")
    if "looseJetCond" not in branchArray:
        rdf = rdf.Define("looseJetCond", "(Jet_pt_JEC > 15.0) && JetIdTight && (Jet_neEmEF + Jet_chEmEF < 0.9)")
        _append_once(branchArray, "looseJetCond")
    if "isLooseJetFailJVM" not in branchArray:
        rdf = rdf.Define("isLooseJetFailJVM", "(looseJetCond && !(PassJetVeto))")
        _append_once(branchArray, "isLooseJetFailJVM")
    if "JVMweight" not in branchArray:
        rdf = rdf.Define("JVMweight", "float(Nonzero(isLooseJetFailJVM).size()==0)")
        _append_once(branchArray, "JVMweight")
    if "Jet_tightNoPt" not in branchArray:
        rdf = rdf.Define("Jet_tightNoPt", "(abs(Jet_eta) < 2.5) && (Jet_rawFactor<0.9) && JetIdTight")
        _append_once(branchArray, "Jet_tightNoPt")
    if "Jet_mediumPtTight" not in branchArray:
        rdf = rdf.Define("Jet_mediumPtTight", "(Jet_pt_JEC > 30.0) && Jet_tightNoPt")
        _append_once(branchArray, "Jet_mediumPtTight")
    if "Jet_drFromMuon" not in branchArray:
        rdf = rdf.Define(
            "Jet_drFromMuon",
            "minDistanceFromEMu(Jet_tightNoPt, Jet_eta, Jet_phi, leadingMuonIdx, leadingElectronIdx, Muon_eta, Muon_phi, Electron_eta, Electron_phi)",
        )
        _append_once(branchArray, "Jet_drFromMuon")
    if "GoodJetCond" not in branchArray:
        rdf = rdf.Define("GoodJetCond", "(Jet_pt_JEC > 30.0) && Jet_tightNoPt && (Jet_drFromMuon>0.4)")
        _append_once(branchArray, "GoodJetCond")

    for jes_source in JES_SOURCES:
        if jes_source not in branchArray:
            continue
        suffix = jes_source.replace("CMS_", "")
        rdf, branchArray = _define_jvm_weight_variation(
            rdf, branchArray, f"{suffix}_up", f"Jet_pt_JEC*(1+{jes_source})"
        )
        rdf, branchArray = _define_jvm_weight_variation(
            rdf, branchArray, f"{suffix}_down", f"Jet_pt_JEC*(1-{jes_source})"
        )

    if all(branch in branchArray for branch in ["Jet_JER_corr", "Jet_JER_corr_up"]):
        rdf, branchArray = _define_jvm_weight_variation(
            rdf, branchArray, "JER_corr_up", "Jet_pt_JEC/Jet_JER_corr*Jet_JER_corr_up"
        )
    if all(branch in branchArray for branch in ["Jet_JER_corr", "Jet_JER_corr_down"]):
        rdf, branchArray = _define_jvm_weight_variation(
            rdf, branchArray, "JER_corr_down", "Jet_pt_JEC/Jet_JER_corr*Jet_JER_corr_down"
        )

    return rdf, recordedModules, branchArray
