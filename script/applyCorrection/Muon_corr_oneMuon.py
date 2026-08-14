import ROOT
import os
import correctionlib


def _has_column(rdf, branch_name, branch_array):
    if branch_name in branch_array:
        return True
    try:
        return bool(rdf.HasColumn(branch_name))
    except Exception:
        return False


def _define_or_redefine(rdf, branch_name, expression, branch_array):
    if _has_column(rdf, branch_name, branch_array):
        rdf = rdf.Redefine(branch_name, expression)
    else:
        rdf = rdf.Define(branch_name, expression)
        branch_array.append(branch_name)
    return rdf


def processing(rdf, recordedModules, branchArray, era, ds=""):
    single_np_eras = {
        "Run3Summer22NanoAODv12",
        "Run3Summer22EENanoAODv12",
        "Run3Summer23NanoAODv12",
        "Run3Summer23BPixNanoAODv12",
    }
    use_private_single_np = era in single_np_eras

    if "Muon_corr_oneMuon.C" not in recordedModules:
        correctionlib.register_pyroot_binding()
        this_dir = os.path.dirname(os.path.abspath(__file__))
        ROOT.gInterpreter.AddIncludePath(this_dir)
        ROOT.gInterpreter.ProcessLine('#include "Muon_corr_oneMuon.C"')
        recordedModules.append("Muon_corr_oneMuon.C")
    this_dir = os.path.dirname(os.path.abspath(__file__))
    module_dir_cpp = this_dir.replace("\\", "/")
    ROOT.gInterpreter.ProcessLine(f'Muon_corr_oneMuon_init("{era}", "{module_dir_cpp}")')

    muon_sf_pt = "Muon_pt"
    for candidate in ("Muon_pt_corr", "Muon_pt_Rcorr", "Muon_pt"):
        if _has_column(rdf, candidate, branchArray):
            muon_sf_pt = candidate
            break

    if _has_column(rdf, "Muon_promptMVA", branchArray):
        rdf = _define_or_redefine(rdf, "Muon_promptMVA_compat", "Muon_promptMVA", branchArray)
    elif _has_column(rdf, "Muon_mvaTTH", branchArray):
        rdf = _define_or_redefine(rdf, "Muon_promptMVA_compat", "Muon_mvaTTH", branchArray)

    rdf = _define_or_redefine(rdf, "Muon_IDscale", f"MuonIDScale_oneMuon(Muon_eta, {muon_sf_pt}, \"nominal\")", branchArray)
    rdf = _define_or_redefine(rdf, "Muon_IDscale_up", f"MuonIDScale_oneMuon(Muon_eta, {muon_sf_pt}, \"systup\")", branchArray)
    rdf = _define_or_redefine(rdf, "Muon_IDscale_down", f"MuonIDScale_oneMuon(Muon_eta, {muon_sf_pt}, \"systdown\")", branchArray)

    if use_private_single_np:
        rdf = _define_or_redefine(rdf, "Muon_Isoscale", "UnitMuonScale_oneMuon(Muon_eta)", branchArray)
        rdf = _define_or_redefine(rdf, "Muon_Isoscale_up", "UnitMuonScale_oneMuon(Muon_eta)", branchArray)
        rdf = _define_or_redefine(rdf, "Muon_Isoscale_down", "UnitMuonScale_oneMuon(Muon_eta)", branchArray)
        rdf = _define_or_redefine(rdf, "Muon_MVAscale", "UnitMuonScale_oneMuon(Muon_eta)", branchArray)
        rdf = _define_or_redefine(rdf, "Muon_MVAscale_up", "UnitMuonScale_oneMuon(Muon_eta)", branchArray)
        rdf = _define_or_redefine(rdf, "Muon_MVAscale_down", "UnitMuonScale_oneMuon(Muon_eta)", branchArray)

        rdf = _define_or_redefine(
            rdf,
            "MuonScale",
            "EventSingleMuonSF_oneMuon(Muon_IDscale, leadingMuonIdx)",
            branchArray,
        )
        rdf = _define_or_redefine(
            rdf,
            "MuonScale_up",
            "EventSingleMuonSF_oneMuon(Muon_IDscale_up, leadingMuonIdx)",
            branchArray,
        )
        rdf = _define_or_redefine(
            rdf,
            "MuonScale_down",
            "EventSingleMuonSF_oneMuon(Muon_IDscale_down, leadingMuonIdx)",
            branchArray,
        )
    else:
        rdf = _define_or_redefine(rdf, "Muon_Isoscale", "UnitMuonScale_oneMuon(Muon_eta)", branchArray)
        rdf = _define_or_redefine(rdf, "Muon_Isoscale_up", "UnitMuonScale_oneMuon(Muon_eta)", branchArray)
        rdf = _define_or_redefine(rdf, "Muon_Isoscale_down", "UnitMuonScale_oneMuon(Muon_eta)", branchArray)
        rdf = _define_or_redefine(rdf, "Muon_MVAscale", f"MuonMVAScale_oneMuon(Muon_eta, {muon_sf_pt}, \"nominal\")", branchArray)
        rdf = _define_or_redefine(rdf, "Muon_MVAscale_up", f"MuonMVAScale_oneMuon(Muon_eta, {muon_sf_pt}, \"systup\")", branchArray)
        rdf = _define_or_redefine(rdf, "Muon_MVAscale_down", f"MuonMVAScale_oneMuon(Muon_eta, {muon_sf_pt}, \"systdown\")", branchArray)

        rdf = _define_or_redefine(
            rdf,
            "MuonScale",
            "EventMuonSF_oneMuon(Muon_IDscale, Muon_Isoscale, Muon_MVAscale, leadingMuonIdx)",
            branchArray,
        )
        rdf = _define_or_redefine(
            rdf,
            "MuonScale_up",
            "MuonScale + std::sqrt("
            "std::pow(EventMuonSF_oneMuon(Muon_IDscale_up, Muon_Isoscale, Muon_MVAscale, leadingMuonIdx) - MuonScale, 2) + "
            "std::pow(EventMuonSF_oneMuon(Muon_IDscale, Muon_Isoscale, Muon_MVAscale_up, leadingMuonIdx) - MuonScale, 2))",
            branchArray,
        )
        rdf = _define_or_redefine(
            rdf,
            "MuonScale_down",
            "std::max(0.0, MuonScale - std::sqrt("
            "std::pow(MuonScale - EventMuonSF_oneMuon(Muon_IDscale_down, Muon_Isoscale, Muon_MVAscale, leadingMuonIdx), 2) + "
            "std::pow(MuonScale - EventMuonSF_oneMuon(Muon_IDscale, Muon_Isoscale, Muon_MVAscale_down, leadingMuonIdx), 2)))",
            branchArray,
        )

    return rdf, recordedModules, branchArray
