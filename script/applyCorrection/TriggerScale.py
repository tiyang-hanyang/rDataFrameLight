import ROOT
import os
import correctionlib
import re


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


def _ordered_unique(items):
    seen = set()
    result = []
    for item in items:
        if not item:
            continue
        if item in seen:
            continue
        seen.add(item)
        result.append(item)
    return result


def _build_trigger_scale_tags(era):
    tags = []
    if not era:
        return tags

    era = str(era).strip()
    tags.append(era)

    if era.endswith(".json"):
        tags.append(os.path.splitext(os.path.basename(era))[0])

    runiii_match = re.match(r"^(RunIII20\d{2}Summer\d+[A-Za-z]*)(NanoAODv\d+\w*)?$", era)
    if runiii_match:
        tags.append(runiii_match.group(1))
        year_match = re.search(r"20\d{2}", runiii_match.group(1))
        if year_match:
            suffix_match = re.search(r"(EE|BPix)$", runiii_match.group(1))
            if suffix_match:
                tags.append(f"Run{year_match.group(0)}{suffix_match.group(1)}")
            tags.append(f"Run{year_match.group(0)}")

    run3_match = re.match(r"^(Run3Summer(\d{2})(EE|BPix)?)(NanoAODv\d+\w*)?$", era)
    if run3_match:
        tags.append(run3_match.group(1))
        if run3_match.group(3):
            tags.append(f"Run20{run3_match.group(2)}{run3_match.group(3)}")
            tags.append(f"Run3Summer{run3_match.group(2)}")
        tags.append(f"Run20{run3_match.group(2)}")

    data_match = re.match(r"^(Run(20\d{2})([A-Z])?)$", era)
    if data_match:
        tags.append(data_match.group(1))
        tags.append(f"Run{data_match.group(2)}")

    year_match = re.search(r"20\d{2}", era)
    if year_match:
        tags.append(f"Run{year_match.group(0)}")

    if "Summer24" in era:
        tags.append("Summer24")

    return _ordered_unique(tags)


def _resolve_trigger_scale_file(this_dir, era):
    search_dirs = [
        os.path.abspath(os.path.join(this_dir, "trigger_sf")),
        os.path.abspath(os.path.join(this_dir, "..", "..", "correction", "custom", "trigger_sf")),
    ]
    era_candidates = _build_trigger_scale_tags(era)

    filenames = []
    for tag in era_candidates:
        filenames.extend(
            [
                f"{tag}_trigger_scale_factor.corr.json",
                f"{tag}_trigger_scale_factor.json",
                f"{tag}_trigger_sf.corr.json",
                f"{tag}_trigger_sf.json",
            ]
        )

    for corr_dir in search_dirs:
        for file_name in filenames:
            candidate = os.path.join(corr_dir, file_name)
            if os.path.exists(candidate):
                return candidate
    return None


def _as_bool(value):
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def _resolve_trigger_scale_source_era(era, options=None):
    options = options or {}

    env_unify = _as_bool(os.environ.get("TRIGGER_SCALE_UNIFY_ERA"))
    opt_unify = _as_bool(options.get("unify_era", False))
    unify_era = env_unify or opt_unify

    if not unify_era:
        return era

    source_era = (
        options.get("unified_era")
        or options.get("source_era")
        or os.environ.get("TRIGGER_SCALE_UNIFIED_ERA")
    )
    if not source_era:
        raise RuntimeError(
            "[TriggerScale] unify_era is enabled, but no unified/source era was provided. "
            "Set moduleOptions['TriggerScale'] = {'unify_era': True, 'unified_era': 'Run2022'} "
            "or export TRIGGER_SCALE_UNIFIED_ERA."
        )
    return source_era


def _detect_trigger_scale_dimension(correction_file, correction_name):
    cset = correctionlib.CorrectionSet.from_file(correction_file)
    corr = cset[correction_name]
    input_names = [item.name for item in corr.inputs]
    if input_names == ["pt", "variation"]:
        return 1
    if input_names == ["leading_pt", "subleading_pt", "variation"]:
        return 2
    raise RuntimeError(
        "[TriggerScale] Unsupported trigger scale correction inputs: "
        f"{input_names}. Expected ['pt', 'variation'] or "
        "['leading_pt', 'subleading_pt', 'variation']."
    )


def processing(rdf, recordedModules, branchArray, era, ds="", options=None):
    this_dir = os.path.dirname(os.path.abspath(__file__))
    source_era = _resolve_trigger_scale_source_era(era, options=options)
    correction_file = _resolve_trigger_scale_file(this_dir, source_era)
    if correction_file is None:
        print(
            f"[TriggerScale] WARNING: no trigger scale correctionlib found for era {era} "
            f"(source era {source_era}); defining unity weights."
        )
        rdf = _define_or_redefine(rdf, "TriggerScale", "1.0f", branchArray)
        rdf = _define_or_redefine(rdf, "TriggerScale_up", "1.0f", branchArray)
        rdf = _define_or_redefine(rdf, "TriggerScale_down", "1.0f", branchArray)
        return rdf, recordedModules, branchArray

    if "TriggerScale.C" not in recordedModules:
        correctionlib.register_pyroot_binding()
        ROOT.gInterpreter.AddIncludePath(this_dir)
        ROOT.gInterpreter.ProcessLine('#include "TriggerScale.C"')
        recordedModules.append("TriggerScale.C")

    correction_name = "leading_muon_trigger_sf"
    correction_dimension = _detect_trigger_scale_dimension(correction_file, correction_name)
    print(
        f"[TriggerScale] INFO: era {era} uses trigger scale source {source_era} "
        f"from file {correction_file}"
    )
    ROOT.gInterpreter.ProcessLine(
        'TriggerScale_init("{}", "{}", {})'.format(
            correction_file.replace("\\", "/"),
            correction_name,
            "true" if correction_dimension == 2 else "false",
        )
    )

    trigger_pt_branch = "Muon_pt"
    for candidate in ("Muon_pt_corr", "Muon_pt_Rcorr", "Muon_pt"):
        if _has_column(rdf, candidate, branchArray):
            trigger_pt_branch = candidate
            break

    scale_expr_base = (
        f"EventTriggerScaleFromRVec({trigger_pt_branch}, leadingMuonIdx, subleadingMuonIdx, "
    )
    rdf = _define_or_redefine(
        rdf,
        "TriggerScale",
        scale_expr_base + "\"nominal\")",
        branchArray,
    )
    rdf = _define_or_redefine(
        rdf,
        "TriggerScale_up",
        scale_expr_base + "\"up\")",
        branchArray,
    )
    rdf = _define_or_redefine(
        rdf,
        "TriggerScale_down",
        scale_expr_base + "\"down\")",
        branchArray,
    )

    return rdf, recordedModules, branchArray
