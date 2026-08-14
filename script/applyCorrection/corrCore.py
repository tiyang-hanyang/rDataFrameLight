import ROOT
import argparse
import importlib.util
import inspect
import os
import json
import sys

from jobDef import GeneralJob


def run_processing_module(procName, procFunc, rdf, recordedModules, branchArray, era, dataset, moduleOptions=None):
    moduleOptions = moduleOptions or {}
    procOptions = moduleOptions.get(procName, {})
    signature = inspect.signature(procFunc)
    if "options" in signature.parameters:
        return procFunc(rdf, recordedModules, branchArray, era, dataset, options=procOptions)
    return procFunc(rdf, recordedModules, branchArray, era, dataset)


def _get_genweight_sum_from_file(file_path):
    file_in = ROOT.TFile(file_path, "read")
    if not file_in or file_in.IsZombie():
        raise OSError(f"cannot open ROOT file: {file_path}")

    try:
        obj = file_in.Get("genWeightSum")
        if obj and obj.InheritsFrom("TH1"):
            return obj.Integral()
        if obj and hasattr(obj, "GetBinContent"):
            return obj.GetBinContent(1)
        if obj and hasattr(obj, "GetVal"):
            return obj.GetVal()

        runs_tree = file_in.Get("Runs")
        if runs_tree and runs_tree.InheritsFrom("TTree"):
            if runs_tree.GetBranch("genEventSumw"):
                total = 0.0
                for entry in runs_tree:
                    total += float(entry.genEventSumw)
                return total
            runs_df = ROOT.ROOT.RDataFrame("Runs", file_path)
            runs_cols = runs_df.GetColumnNames()
            if "genEventSumw" in runs_cols:
                return runs_df.Sum("genEventSumw").GetValue()

        events_tree = file_in.Get("Events")
        if events_tree and events_tree.InheritsFrom("TTree"):
            if events_tree.GetBranch("genWeight"):
                total = 0.0
                for entry in events_tree:
                    total += float(entry.genWeight)
                return total
            events_df = ROOT.ROOT.RDataFrame("Events", file_path)
            event_cols = events_df.GetColumnNames()
            if "genWeight" in event_cols:
                return events_df.Sum("genWeight").GetValue()

        key_names = [key.GetName() for key in file_in.GetListOfKeys()] if file_in.GetListOfKeys() else []
        runs_branches = []
        if runs_tree and runs_tree.InheritsFrom("TTree"):
            runs_branches = [br.GetName() for br in runs_tree.GetListOfBranches()]
        event_branches = []
        if events_tree and events_tree.InheritsFrom("TTree"):
            event_branches = [br.GetName() for br in events_tree.GetListOfBranches()]

        raise RuntimeError(
            f"{file_path} has no usable genWeight source; "
            f"keys={key_names}, Runs branches={runs_branches[:10]}, Events has genWeight={'genWeight' in event_branches}"
        )
    finally:
        file_in.Close()


def _prepare_snapshot_columns(branchArray):
    if branchArray is None:
        raise RuntimeError("branchArray is None before Snapshot")

    normalized = []
    bad_entries = []
    for idx, branch in enumerate(branchArray):
        if branch is None:
            bad_entries.append((idx, None, "NoneType"))
            continue
        try:
            branch_name = str(branch)
        except Exception as exc:
            bad_entries.append((idx, repr(branch), f"{type(branch).__name__}: {exc}"))
            continue
        if not branch_name:
            bad_entries.append((idx, repr(branch), "empty string"))
            continue
        normalized.append(branch_name)

    if bad_entries:
        raise RuntimeError(f"Invalid Snapshot branch entries: {bad_entries[:10]}")
    if not normalized:
        raise RuntimeError("No branches selected for Snapshot")

    out = ROOT.std.vector("string")()
    for branch_name in normalized:
        out.push_back(branch_name)
    return out


def _chunk_file_paths(file_paths, chunk_size):
    if chunk_size is None:
        return [file_paths]
    chunk_size = int(chunk_size)
    if chunk_size <= 0:
        raise ValueError(f"heavySplit chunk size must be positive, got {chunk_size}")
    return [file_paths[i:i + chunk_size] for i in range(0, len(file_paths), chunk_size)]

neededBr = [
    "run",
    "luminosityBlock",
    "event",

    "PV_npvsGood",
    "Pileup_nPU",
    "Pileup_nTrueInt",

    "LHEPdfWeight",
    "LHEScaleWeight",
    "PSWeight",
    "genWeight",

    "nGenPart",
    "GenPart_genPartIdxMother",
    "GenPart_status",
    "GenPart_statusFlags",
    "GenPart_pdgId",
    "GenPart_pt",
    "GenPart_eta",
    "GenPart_phi",
    "GenPart_mass",

    "nJet",
    "Jet_jetId",
    "Jet_area",
    "Jet_pt",
    "Jet_eta",
    "Jet_phi",
    "Jet_mass",
    "Jet_btagPNetB",
    "Jet_btagRobustParTAK4B",
    "Jet_btagUParTAK4B",
    "Jet_neHEF",
    "Jet_neEmEF",
    "Jet_nConstituents",
    "Jet_muEF",
    "Jet_chHEF",
    "Jet_chMultiplicity",
    "Jet_chEmEF",
    "Jet_neMultiplicity",
    "Jet_hadronFlavour",

    "nGoodJet",
    "GoodJetCond",
    "leadingJetIdx",
    "subleadingJetIdx",
    "thirdJetIdx",
    "nBJet",
    "BJetCond",
    "leadingBJetIdx",
    "subleadingBJetIdx",
    "thirdBJetIdx",

    "nMuon",
    "Muon_pt",
    "Muon_ptErr",
    "Muon_eta",
    "Muon_phi",
    "Muon_mass",
    "Muon_charge",
    "Muon_miniPFRelIso_all",
    "Muon_pfRelIso04_all",
    "Muon_pfIsoId",
    "Muon_miniIsoId",
    "Muon_looseId",
    "Muon_mediumId",
    "Muon_tightId",
    "Muon_mvaMuID_WP",
    "Muon_mvaMuID",
    "Muon_promptMVA",
    "Muon_promptMVA_compat",
    "Muon_mvaTTH",
    "Muon_dxy",
    "Muon_dxyErr",
    "Muon_dz",
    "Muon_dzErr",
    "Muon_sip3d",
    "Muon_jetDF",
    "Muon_jetIdx",
    "Muon_jetPtRelv2",
    "Muon_jetRelIso",
    "Muon_nTrackerLayers",

    "isGoodMuon_mva",
    "nGoodMuon_mva",
    "leadingMuonIdx",
    "subleadingMuonIdx",

    "Muon_pt_Rscale",
    "Muon_pt_Rcorr",
    "Muon_pt_Rscale_up",
    "Muon_pt_Rscale_dn",
    "Muon_pt_Rcorr_resolup",
    "Muon_pt_Rcorr_resoldn",
    "Muon_originTraceCode",
    "Muon_hasGenMuon",
    "Muon_isFromB",
    "Muon_isFromC",
    "Muon_isPrompt",
    "Muon_isFromFromTop",

    "leadingMuon_charge",
    "leadingMuon_pt",
    "leadingMuon_eta", 
    "leadingMuon_phi", 
    "leadingMuon_mass", 
    "leadingMuonP4",
    "subleadingMuon_charge", 
    "subleadingMuon_pt", 
    "subleadingMuon_eta",
    "subleadingMuon_phi", 
    "subleadingMuon_mass", 
    "subleadingMuonP4",
    "diMuonP4", 
    "diMuon_mass", 
    "diMuon_deltaR", 
    "diMuon_deltaPhi", 
    "diMuon_deltaEta", 

    "MET_phi",
    "MET_pt",
    "MET_sumEt",
    "PuppiMET_phi",
    "PuppiMET_pt",
    "PuppiMET_sumEt",

    "HLT_Mu3_PFJet40",
    "HLT_Mu8",
    "HLT_Mu17",
    "HLT_Mu20",
    "HLT_Mu27", 
    "HLT_IsoMu20",
    "HLT_IsoMu24",
    "HLT_IsoMu24_eta2p1",
    "HLT_IsoMu27",
    "HLT_Mu50",
    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
    "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8",
    "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8",
    "HLT_TripleMu_5_3_3_Mass3p8_DZ",
    "HLT_TripleMu_10_5_5_DZ",
    "HLT_TripleMu_12_10_5",

    "Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    "Flag_eeBadScFilter",
    "Flag_BadPFMuonDzFilter",
    "Flag_hfNoisyHitsFilter",

    "HLT_CascadeMu100",
    "HLT_HighPtTkMu100",

    "weight_XS",
    "PassJetVeto",
    "btag_weight",
    "btag_weight_up",
    "btag_weight_down",
    "btag_weight_lf_up",
    "btag_weight_lf_down",
    "btag_weight_lfstats1_up",
    "btag_weight_lfstats1_down",
    "btag_weight_lfstats2_up",
    "btag_weight_lfstats2_down",
    "btag_weight_hf_up",
    "btag_weight_hf_down",
    "btag_weight_hfstats1_up",
    "btag_weight_hfstats1_down",
    "btag_weight_hfstats2_up",
    "btag_weight_hfstats2_down",
    "btag_weight_cferr1_up",
    "btag_weight_cferr1_down",
    "btag_weight_cferr2_up",
    "btag_weight_cferr2_down",
    "btag_weight_mur_up",
    "btag_weight_mur_down",
    "btag_weight_pileup_up",
    "btag_weight_pileup_down",
    "btag_weight_pdfas_up",
    "btag_weight_pdfas_down",
    "btag_weight_isrdef_up",
    "btag_weight_isrdef_down",
    "btag_weight_topmass_up",
    "btag_weight_topmass_down",
    "btag_weight_bfragmentation_up",
    "btag_weight_bfragmentation_down",
    "btag_weight_hdamp_up",
    "btag_weight_hdamp_down",
    "btag_weight_jer_up",
    "btag_weight_jer_down",
    "btag_weight_statistic_up",
    "btag_weight_statistic_down",
    "btag_weight_muf_up",
    "btag_weight_muf_down",
    "btag_weight_fsrdef_up",
    "btag_weight_fsrdef_down",
    "btag_weight_jes_up",
    "btag_weight_jes_down",
    "btag_weight_type3_up",
    "btag_weight_type3_down",
    "Muon_IDscale",
    "Muon_IDscale_up",
    "Muon_IDscale_down",
    "Muon_Isoscale",
    "Muon_Isoscale_up",
    "Muon_Isoscale_down",
    "Muon_MVAscale",
    "Muon_MVAscale_up",
    "Muon_MVAscale_down",
    "MuonScale",
    "MuonScale_up",
    "MuonScale_down",
    "TriggerScale",
    "TriggerScale_up",
    "TriggerScale_down",
    "ZptWgt",
    "PUWeight",
    "PUWeight_up",
    "PUWeight_down",
    "IsoMu24Scale",

    "Jet_rawFactor",
    "Jet_JER_corr",
    "Jet_JER_corr_up",
    "Jet_JER_corr_down",
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
    "Rho_fixedGridRhoFastjetAll",
    "Jet_genJetIdx",
    "nGenJet",
    "GenJet_pt",
    "GenJet_eta",
    "GenJet_phi",
    "GenJet_hadronFlavour",
    "RawPuppiMET_phi",
    "RawPuppiMET_pt",
    "RawPuppiMET_sumEt",
    "Jet_muonSubtrFactor",

    "nonPromptWeight",
    "Muon_conePt",
    "onlyLeadingFake",
    "onlySubleadingFake",
    "bothFake",

    "Electron_pt",
    "Electron_phi",
    "Electron_eta",
    "Electron_superclusterEta",
    "Electron_deltaEtaSC",
    "Electron_sieie",
    "Electron_mass",
    "Electron_charge",

    "Electron_genPartIdx",
    "Electron_genPartFlav",

    "Electron_dxy",
    "Electron_dz",
    "Electron_ip3d",
    "Electron_sip3d",

    "Electron_eInvMinusPInv",
    "Electron_hoe",
    "Electron_convVeto",
    "Electron_lostHits",
    "Electron_mvaIso_WP80",
    "Electron_mvaIso_WP90",
    "Electron_mvaIso_WPHZZ",
    "Electron_mvaNoIso_WP80",
    "Electron_mvaNoIso_WP90",
    "Electron_jetIdx",
    "Electron_jetRelIso",
    "Electron_jetPtRelv2",
    "Electron_jetDF",
    "Electron_miniPFRelIso_all",
    "Electron_pfRelIso03_all",
    "Electron_pfRelIso04_all",
    "Electron_cutBased",
    "Electron_promptMVA",

    "leadingElectronIdx"
]

_jvm_syst_sources = [
    "scale_j_FlavorQCD",
    "scale_j_RelativeBal",
    "scale_j_HF",
    "scale_j_BBEC1",
    "scale_j_EC2",
    "scale_j_Absolute",
    "scale_j_Absolute_2022",
    "scale_j_HF_2022",
    "scale_j_EC2_2022",
    "scale_j_RelativeSample_2022",
    "scale_j_BBEC1_2022",
    "scale_j_Absolute_2022EE",
    "scale_j_HF_2022EE",
    "scale_j_EC2_2022EE",
    "scale_j_RelativeSample_2022EE",
    "scale_j_BBEC1_2022EE",
    "scale_j_Absolute_2023",
    "scale_j_HF_2023",
    "scale_j_EC2_2023",
    "scale_j_RelativeSample_2023",
    "scale_j_BBEC1_2023",
    "scale_j_Absolute_2023BPix",
    "scale_j_HF_2023BPix",
    "scale_j_EC2_2023BPix",
    "scale_j_RelativeSample_2023BPix",
    "scale_j_BBEC1_2023BPix",
    "scale_j_Absolute_2024",
    "scale_j_HF_2024",
    "scale_j_EC2_2024",
    "scale_j_RelativeSample_2024",
    "scale_j_BBEC1_2024",
]
for _jvm_syst_source in _jvm_syst_sources:
    neededBr.append(f"JVMweight_{_jvm_syst_source}_up")
    neededBr.append(f"JVMweight_{_jvm_syst_source}_down")
neededBr.append("JVMweight_JER_corr_up")
neededBr.append("JVMweight_JER_corr_down")

# to process all the decorations for merged process into a single output
def processMergeDS(
    era,
    dataset,
    filePaths,
    commonOutDir,
    procedures,
    recordedModules,
    needSlice=1,
    moduleOptions=None,
    output_name=None,
):
    outDir = commonOutDir+"/"+era+"/"+dataset+"/"
    if not os.path.exists(outDir):
        os.makedirs(outDir)
        print("create out folder")
    # always enable MT
    ROOT.EnableImplicitMT() 

    isData = (era=="Run2024C") or (era=="Run2024D") or (era=="Run2024E") or (era=="Run2024F") or (era=="Run2024G") or (era=="Run2024H") or (era=="Run2024I") or (era=="Run2023C") or (era=="Run2023D") or (era=="Run2022C") or (era=="Run2022D") or (era=="Run2022E") or (era=="Run2022F") or (era=="Run2022G")

    # if the merged job output already exist, then nothing needed, just skip.
    if output_name is None:
        output_name = dataset + "_skimmed.root"
    fout = outDir + output_name
    if os.path.isfile(fout):
        print(fout, "already exists, skip")
        return recordedModules

    # processing
    ch1=ROOT.TChain("Events")
    samples = 0
    validGenWeightSums = []
    for fin in filePaths:
        fileTest = ROOT.TFile(fin, "read")
        if not "Events" in fileTest.GetListOfKeys():
            fileTest.Close()
            continue
        tempTree = fileTest.Get("Events")
        if not tempTree:
            fileTest.Close()
            continue
        if tempTree.GetEntries() == 0:
            fileTest.Close()
            continue
        fileTest.Close()

        if not isData:
            fileGenWeightSum = _get_genweight_sum_from_file(fin)
            validGenWeightSums.append(fileGenWeightSum)

        ch1.Add(fin)
        samples += 1
    if samples == 0:
        return recordedModules
    rdf = ROOT.ROOT.RDataFrame(ch1)
    opt = ROOT.RDF.RSnapshotOptions()
    if needSlice:
        branchArray = list([i for i in neededBr if i in rdf.GetColumnNames()])  
    else:
        branchArray = list([i for i in rdf.GetColumnNames()])  
    this_dir = os.path.dirname(os.path.abspath(__file__))
    for (procName, procedureFile) in procedures:
        module_name = f"processing_{procName}"
        module_path = this_dir+"/"+procedureFile
        spec = importlib.util.spec_from_file_location(module_name, module_path)
        if spec is None or spec.loader is None:
            raise ImportError(f"Cannot load spec from: {module_path}")
        procModel = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(procModel)
        if not hasattr(procModel, "processing"):
            raise AttributeError(f"{module_path} does not define function `processing`")
        procFunc = getattr(procModel, "processing")
        if not callable(procFunc):
            raise TypeError(f"`processing` in {module_path} is not callable")
        rdf, recordedModules, branchArray = run_processing_module(
            procName,
            procFunc,
            rdf,
            recordedModules,
            branchArray,
            era,
            dataset,
            moduleOptions=moduleOptions,
        )
    snapshotColumns = _prepare_snapshot_columns(branchArray)
    rdf.Snapshot("Events", fout, snapshotColumns, opt)

    # transfer the genWeight
    # only for MC
    if not isData:
        totalGenWeight = sum(validGenWeightSums)
        fileOut = ROOT.TFile(fout, "UPDATE")
        genWeightSumHist = ROOT.TH1D("genWeightSum", "sum of genWeight", 1, 0.0, 1.0) 
        genWeightSumHist.SetBinContent(1, totalGenWeight)
        genWeightSumHist.GetYaxis().SetTitle("sum(genWeight)")
        fileOut.Write("")
        fileOut.Close()

    return recordedModules

# to process all the decorations for the non-merged process, one file to one file
def processNonMergeDS(era, dataset, filePaths, commonOutDir, procedures, recordedModules, batch_size=16, needSlice=1, moduleOptions=None):
    outDir = commonOutDir+"/"+era+"/"+dataset+"/"
    if not os.path.exists(outDir):
        os.makedirs(outDir)
        print("create out folder")
    # always enable MT
    ROOT.EnableImplicitMT() 

    # define jobs IO
    jobs = []
    for fin in filePaths:
        fout = outDir + fin.split("/")[-1].split(".root")[0] + "_skimmed.root"
        if not os.path.isfile(fout):
            jobs.append((fin, fout))
        else:
            print(fout, "already eixsts")
    if not jobs:
        print("All outputs already exist.")
        return recordedModules
    print(len(jobs), "files to process")
    
    # shared option
    opt = ROOT.RDF.RSnapshotOptions()
    opt.fLazy = True

    # creat batch run
    print("batch_size:", batch_size)
    for s in range(0, len(jobs), batch_size):
        sub = jobs[s:s+batch_size]
        graphs, handles = [], []
        for fin, fout in sub:
            rdf = ROOT.ROOT.RDataFrame("Events", fin)
            if needSlice:
                branchArray = list([i for i in neededBr if i in rdf.GetColumnNames()])  
            else:
                branchArray = list([i for i in rdf.GetColumnNames()])  
            this_dir = os.path.dirname(os.path.abspath(__file__))
            for (procName, procedureFile) in procedures:
                module_name = f"processing_{procName}"
                module_path = this_dir+"/"+procedureFile
                spec = importlib.util.spec_from_file_location(module_name, module_path)
                if spec is None or spec.loader is None:
                    raise ImportError(f"Cannot load spec from: {module_path}")
                procModel = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(procModel)
                if not hasattr(procModel, "processing"):
                    raise AttributeError(f"{module_path} does not define function `processing`")
                procFunc = getattr(procModel, "processing")
                if not callable(procFunc):
                    raise TypeError(f"`processing` in {module_path} is not callable")
                rdf, recordedModules, branchArray = run_processing_module(
                    procName,
                    procFunc,
                    rdf,
                    recordedModules,
                    branchArray,
                    era,
                    dataset,
                    moduleOptions=moduleOptions,
                )
            snapshotColumns = _prepare_snapshot_columns(branchArray)
            h = rdf.Snapshot("Events", fout, snapshotColumns, opt)
            handles.append(h)  
            graphs.append(rdf)
        ROOT.RDF.RunGraphs(handles)
        for h in handles:
            h.GetValue() 
        print(f"batch collected: {s+len(sub)}/{len(jobs)}")
        for fin, fout in sub:
            fileIn = ROOT.TFile(fileIn, "read")
            genWeightHist = fileIn.Get("genWeightSum")
            genWeightHist.SetDirectory(0)
            fileIn.Close()
            fileOut = ROOT.TFile(fileOut, "UPDATE")
            fileOut.Write("")
            fileOut.Close()
        print(f"batch done: {s+len(sub)}/{len(jobs)}")
    return recordedModules

# get the job configuration without pre-registration
def load_job_from_file(file_path: str) -> GeneralJob:
    file_path = os.path.abspath(file_path)
    if not os.path.isfile(file_path):
        raise FileNotFoundError(file_path)
    module_name = "_dynamic_job_module"
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    for name, obj in vars(module).items():
        if inspect.isclass(obj) and issubclass(obj, GeneralJob) and obj is not GeneralJob:
            return obj() 
    raise RuntimeError(f"{file_path} is not a derived class of GeneralJob")


def _load_sample_list(job, era, ds):
    periods, datasets, mergeDS, workflow, fileJson, outDir = job.declare()
    if era not in periods:
        raise RuntimeError(f"Era {era} is not declared in job periods: {periods}")
    if ds not in datasets:
        raise RuntimeError(f"Dataset {ds} is not declared in job datasets: {datasets}")

    with open(fileJson[era]) as jFile:
        jsonFull = json.load(jFile)

    if ds not in jsonFull["dir"].keys() or ds not in jsonFull["file"].keys():
        raise RuntimeError(f"{ds} is not available in sample json for era {era}")

    sampleList = sorted([jsonFull["dir"][ds] + i for i in jsonFull["file"][ds]])
    if len(sampleList) == 0:
        raise RuntimeError(f"Dataset {ds} in era {era} has an empty sample list")
    return sampleList


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "job_file",
        help="write a job file to inherit the jobDef.py and give here "
    )
    args = parser.parse_args()
    job = load_job_from_file(args.job_file)
    periods, datasets, mergeDS, workflow, fileJson, outDir = job.declare()
    moduleOptions = getattr(job, "moduleOptions", {})
    heavySplit = getattr(job, "heavySplit", {})

    # need additional control only turned on the reduce branches in the initial skim
    needSlice = 0

    # to avoid load already loaded macros 
    recordedModules=[]
    for era in periods:
        print("skim era: " + str(era))
        with open(fileJson[era]) as jFile:
            jsonFull = json.load(jFile)
        for ds in datasets: 
            print("skim ds: "+ str(ds))
            if ds not in jsonFull["dir"].keys() or ds not in jsonFull["file"].keys():
                print(ds, "not available after selection.")
                continue
            sampleList = sorted([jsonFull["dir"][ds] + i for i in jsonFull["file"][ds]])
            if len(sampleList) == 0:
                print("DS", ds, "does not exist")
                continue
            print(sampleList)
            if ds in mergeDS:
                chunk_size = heavySplit.get(ds)
                if chunk_size is None:
                    recordedModules = processMergeDS(
                        era,
                        ds,
                        sampleList,
                        outDir,
                        workflow,
                        recordedModules,
                        needSlice,
                        moduleOptions=moduleOptions,
                    )
                else:
                    sample_chunks = _chunk_file_paths(sampleList, chunk_size)
                    print(f"heavySplit enabled for {ds}: {len(sampleList)} input files -> {len(sample_chunks)} merged chunks with size {chunk_size}")
                    for chunk_idx, sample_chunk in enumerate(sample_chunks, start=1):
                        output_name = f"{ds}_skimmed_{int(chunk_size)}_{chunk_idx}.root"
                        print(f"processing heavySplit chunk {chunk_idx}/{len(sample_chunks)}: {output_name}")
                        recordedModules = processMergeDS(
                            era,
                            ds,
                            sample_chunk,
                            outDir,
                            workflow,
                            recordedModules,
                            needSlice,
                            moduleOptions=moduleOptions,
                            output_name=output_name,
                        )
            else:
                recordedModules=processNonMergeDS(era, ds, sampleList, outDir, workflow, recordedModules, needSlice, moduleOptions=moduleOptions)

if __name__ == "__main__":
    main()
