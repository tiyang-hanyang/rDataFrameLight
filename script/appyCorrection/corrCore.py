import ROOT
import argparse
import importlib.util
import inspect
import os
import json

from jobDef import GeneralJob

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
    "Muon_IDscale",
    "Muon_Isoscale",
    "Muon_MVAscale",
    "MuonScale",
    "ZptWgt",
    "PUWeight",
    "IsoMu24Scale",

    "Jet_rawFactor",
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
    "bothFake"
]

# to process all the decorations for merged process into a single output
def processMergeDS(era, dataset, filePaths, commonOutDir, procedures, recordedModules, needSlice=1):
    outDir = commonOutDir+"/"+era+"/"+dataset+"/"
    if not os.path.exists(outDir):
        os.makedirs(outDir)
        print("create out folder")
    # always enable MT
    ROOT.EnableImplicitMT() 

    # if the merged job output already exist, then nothing needed, just skip.
    fout = outDir + dataset + "_skimmed.root"
    if os.path.isfile(fout):
        print(fout, "already exists, skip")
        return recordedModules

    # processing
    ch1=ROOT.TChain("Events")
    for fin in filePaths:
        fileTest = ROOT.TFile(fin, "read")
        if not "Events" in fileTest.GetListOfKeys():
            continue
        ch1.Add(fin)
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
        rdf, recordedModules, branchArray = procFunc(rdf, recordedModules, branchArray, era, dataset)
    rdf.Snapshot("Events", fout, branchArray, opt)

    # transfer the genWeight
    # only for MC
    isData = (era=="Run2024C") or (era=="Run2024D") or (era=="Run2024E") or (era=="Run2024F") or (era=="Run2024G") or (era=="Run2024H") or (era=="Run2024I") or (era=="Run2023C") or (era=="Run2023D") or (era=="Run2022C") or (era=="Run2022D") or (era=="Run2022E") or (era=="Run2022F") or (era=="Run2022G")
    if not isData:
        totalGenWeight = 0.0
        for fin in filePaths:
            fileIn = ROOT.TFile(fin, "read")
            genWeightHist = fileIn.Get("genWeightSum")
            totalGenWeight += genWeightHist.Integral()
            fileIn.Close()
        fileOut = ROOT.TFile(fout, "UPDATE")
        genWeightSumHist = ROOT.TH1D("genWeightSum", "sum of genWeight", 1, 0.0, 1.0) 
        genWeightSumHist.SetBinContent(1, totalGenWeight)
        genWeightSumHist.GetYaxis().SetTitle("sum(genWeight)")
        fileOut.Write("")
        fileOut.Close()

    return recordedModules

# to process all the decorations for the non-merged process, one file to one file
def processNonMergeDS(era, dataset, filePaths, commonOutDir, procedures, recordedModules, batch_size=16, needSlice=1):
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
                rdf, recordedModules, branchArray = procFunc(rdf, recordedModules, branchArray, era, dataset)
            h = rdf.Snapshot("Events", fout, branchArray, opt)
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "job_file",
        help="write a job file to inherit the jobDef.py and give here "
    )
    args = parser.parse_args()
    job = load_job_from_file(args.job_file)
    periods, datasets, mergeDS, workflow, fileJson, outDir = job.declare()

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
            sampleList = [ jsonFull["dir"][ds] + i for i in jsonFull["file"][ds]]
            print(sampleList)
            if ds in mergeDS:
                recordedModules=processMergeDS(era, ds, sampleList, outDir, workflow, recordedModules, needSlice)
            else:
                recordedModules=processNonMergeDS(era, ds, sampleList, outDir, workflow, recordedModules, needSlice)

if __name__ == "__main__":
    main()
