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

# muon scale correction only now, later consider about the corrections for the trigger
def processing(rdf, recordedModules, branchArray, era, ds=""):
    single_np_eras = {
        "Run3Summer22NanoAODv12",
        "Run3Summer22EENanoAODv12",
        "Run3Summer23NanoAODv12",
        "Run3Summer23BPixNanoAODv12",
    }
    use_private_single_np = era in single_np_eras

    if "Muon_corr.C" not in recordedModules:
        correctionlib.register_pyroot_binding()
        this_dir = os.path.dirname(os.path.abspath(__file__))
        ROOT.gInterpreter.AddIncludePath(this_dir)
        ROOT.gInterpreter.ProcessLine('#include "Muon_corr.C"')
        recordedModules.append("Muon_corr.C")
    this_dir = os.path.dirname(os.path.abspath(__file__))
    module_dir_cpp = this_dir.replace("\\", "/")
    ROOT.gInterpreter.ProcessLine(f'Muon_corr_init("{era}", "{module_dir_cpp}")')

    # Evaluate muon SFs on the best available corrected muon pT without
    # overwriting the original NanoAOD Muon_pt branch.
    muon_sf_pt = "Muon_pt"
    for candidate in ("Muon_pt_corr", "Muon_pt_Rcorr", "Muon_pt"):
        if _has_column(rdf, candidate, branchArray):
            muon_sf_pt = candidate
            break

    # v15 stores the prompt-MVA branch as Muon_promptMVA, while v12 uses
    # Muon_mvaTTH for the same role. Provide one stable compatibility name for
    # downstream cut/plot configs.
    if _has_column(rdf, "Muon_promptMVA", branchArray):
        rdf = _define_or_redefine(rdf, "Muon_promptMVA_compat", "Muon_promptMVA", branchArray)
    elif _has_column(rdf, "Muon_mvaTTH", branchArray):
        rdf = _define_or_redefine(rdf, "Muon_promptMVA_compat", "Muon_mvaTTH", branchArray)

    if use_private_single_np:
        rdf = _define_or_redefine(rdf, "Muon_IDscale", f"MuonIDScale(Muon_eta, {muon_sf_pt}, \"nominal\")", branchArray)
        rdf = _define_or_redefine(rdf, "Muon_IDscale_up", f"MuonIDScale(Muon_eta, {muon_sf_pt}, \"systup\")", branchArray)
        rdf = _define_or_redefine(rdf, "Muon_IDscale_down", f"MuonIDScale(Muon_eta, {muon_sf_pt}, \"systdown\")", branchArray)

        # Keep the legacy component branches available while collapsing the fit to
        # a single muon-scale nuisance parameter.
        rdf = _define_or_redefine(rdf, "Muon_Isoscale", "UnitMuonScale(Muon_eta)", branchArray)
        rdf = _define_or_redefine(rdf, "Muon_Isoscale_up", "UnitMuonScale(Muon_eta)", branchArray)
        rdf = _define_or_redefine(rdf, "Muon_Isoscale_down", "UnitMuonScale(Muon_eta)", branchArray)
        rdf = _define_or_redefine(rdf, "Muon_MVAscale", "UnitMuonScale(Muon_eta)", branchArray)
        rdf = _define_or_redefine(rdf, "Muon_MVAscale_up", "UnitMuonScale(Muon_eta)", branchArray)
        rdf = _define_or_redefine(rdf, "Muon_MVAscale_down", "UnitMuonScale(Muon_eta)", branchArray)

        rdf = _define_or_redefine(
            rdf,
            "MuonScale",
            f"EventSingleMuonSFSelected(Muon_IDscale, isGoodMuon_mva, {muon_sf_pt})",
            branchArray,
        )
        rdf = _define_or_redefine(
            rdf,
            "MuonScale_up",
            f"EventSingleMuonSFSelected(Muon_IDscale_up, isGoodMuon_mva, {muon_sf_pt})",
            branchArray,
        )
        rdf = _define_or_redefine(
            rdf,
            "MuonScale_down",
            f"EventSingleMuonSFSelected(Muon_IDscale_down, isGoodMuon_mva, {muon_sf_pt})",
            branchArray,
        )
    else:
        rdf = _define_or_redefine(rdf, "Muon_IDscale", f"MuonIDScale(Muon_eta, {muon_sf_pt}, \"nominal\")", branchArray)
        rdf = _define_or_redefine(rdf, "Muon_IDscale_up", f"MuonIDScale(Muon_eta, {muon_sf_pt}, \"systup\")", branchArray)
        rdf = _define_or_redefine(rdf, "Muon_IDscale_down", f"MuonIDScale(Muon_eta, {muon_sf_pt}, \"systdown\")", branchArray)

        # 2024 keeps the official ID and prompt-MVA corrections, but excludes the
        # loose mini-iso term entirely. The fit-facing muon-scale NP is therefore
        # built from ID + WP64 MVA only, with no loose-iso contribution in either
        # the nominal MuonScale or its up/down variations.
        rdf = _define_or_redefine(rdf, "Muon_Isoscale", "UnitMuonScale(Muon_eta)", branchArray)
        rdf = _define_or_redefine(rdf, "Muon_Isoscale_up", "UnitMuonScale(Muon_eta)", branchArray)
        rdf = _define_or_redefine(rdf, "Muon_Isoscale_down", "UnitMuonScale(Muon_eta)", branchArray)

        rdf = _define_or_redefine(rdf, "Muon_MVAscale", f"MuonMVAScale(Muon_eta, {muon_sf_pt}, \"nominal\")", branchArray)
        rdf = _define_or_redefine(rdf, "Muon_MVAscale_up", f"MuonMVAScale(Muon_eta, {muon_sf_pt}, \"systup\")", branchArray)
        rdf = _define_or_redefine(rdf, "Muon_MVAscale_down", f"MuonMVAScale(Muon_eta, {muon_sf_pt}, \"systdown\")", branchArray)

        rdf = _define_or_redefine(
            rdf,
            "MuonScale",
            f"EventMuonSFSelected(Muon_IDscale, Muon_Isoscale, Muon_MVAscale, isGoodMuon_mva, {muon_sf_pt})",
            branchArray,
        )
        rdf = _define_or_redefine(
            rdf,
            "MuonScale_up",
            "MuonScale + std::sqrt("
            f"std::pow(EventMuonSFSelected(Muon_IDscale_up, Muon_Isoscale, Muon_MVAscale, isGoodMuon_mva, {muon_sf_pt}) - MuonScale, 2) + "
            f"std::pow(EventMuonSFSelected(Muon_IDscale, Muon_Isoscale, Muon_MVAscale_up, isGoodMuon_mva, {muon_sf_pt}) - MuonScale, 2))",
            branchArray,
        )
        rdf = _define_or_redefine(
            rdf,
            "MuonScale_down",
            "std::max(0.0, MuonScale - std::sqrt("
            f"std::pow(MuonScale - EventMuonSFSelected(Muon_IDscale_down, Muon_Isoscale, Muon_MVAscale, isGoodMuon_mva, {muon_sf_pt}), 2) + "
            f"std::pow(MuonScale - EventMuonSFSelected(Muon_IDscale, Muon_Isoscale, Muon_MVAscale_down, isGoodMuon_mva, {muon_sf_pt}), 2)))",
            branchArray,
        )

    return rdf, recordedModules, branchArray
