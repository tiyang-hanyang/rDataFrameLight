import ROOT
import os

# redefine the branch operation method and directory to save

def split(era, dataset, filePaths):
    # first check the existence of output!
    outDir = "splitPrompt/"+era+"/"+dataset+"/"
    if not os.path.exists(outDir):
        os.makedirs(outDir)
        print("create out folder")

    outLeadingPromptFile = "splitPrompt/"+era+"/"+dataset+"/"+dataset+"_leadingMuonPrompt_incluTau.root"
    outLeadingNonFile = "splitPrompt/"+era+"/"+dataset+"/"+dataset+"_leadingMuonNonPrompt_incluTau.root"
    outSubLeadingPromptFile = "splitPrompt/"+era+"/"+dataset+"/"+dataset+"_subleadingMuonPrompt_incluTau.root"
    outSubLeadingNonFile = "splitPrompt/"+era+"/"+dataset+"/"+dataset+"_subleadingMuonNonPrompt_incluTau.root"

    skipLeading = 0
    skipSubleading = 0
    if os.path.isfile(outLeadingNonFile) and os.path.isfile(outLeadingPromptFile):
        skipLeading = 1
    if os.path.isfile(outSubLeadingNonFile) and os.path.isfile(outSubLeadingPromptFile):
        skipSubleading = 1

    if skipLeading and skipSubleading:
        return
    
    print("spliting: ", dataset)

    ROOT.EnableImplicitMT()

    # run merged as not large
    rdf = ROOT.RDataFrame("Events", filePaths)

    rdf = rdf.Define("GenMuonCond", "abs(GenPart_pdgId)==13 && GenPart_status==1")
    rdf = rdf.Define("GenMuonIdx", "Nonzero(GenMuonCond)")
    # rdf = rdf.Define("PromptMuonCond", "abs(GenPart_pdgId)==13 && GenPart_status==1 && (GenPart_statusFlags & (1<<0))")
    # add the tau definition inside
    rdf = rdf.Define("PromptMuonCond", "abs(GenPart_pdgId)==13 && GenPart_status==1 && ( (GenPart_statusFlags & (1<<0)) || (GenPart_statusFlags & (1<<3)) )")
    rdf = rdf.Define("leadingMuonEta", "Muon_eta[leadingMuonIdx]").Define("leadingMuonPhi", "Muon_phi[leadingMuonIdx]").Define("subleadingMuonEta", "Muon_eta[subleadingMuonIdx]").Define("subleadingMuonPhi", "Muon_phi[subleadingMuonIdx]")
    rdf = rdf.Define("leadingGenMuonIdx",
        "int bestIdx = -1; \
        float bestDR = 0.2f;  \
        for (auto i : GenMuonIdx) { \
        float dr = ROOT::VecOps::DeltaR(leadingMuonEta, GenPart_eta[i], leadingMuonPhi, GenPart_phi[i]);    \
        if (dr < bestDR) {bestDR = dr; bestIdx = (int)i;}} \
        return bestIdx;"
    )
    rdf = rdf.Define("subleadingGenMuonIdx",
        "int bestIdx = -1; \
        float bestDR = 0.2f; \
        for (auto i : GenMuonIdx) { \
            float dr = ROOT::VecOps::DeltaR(subleadingMuonEta, GenPart_eta[i], subleadingMuonPhi, GenPart_phi[i]);    \
            if (dr < bestDR) {bestDR = dr; bestIdx = (int)i;}} \
        return bestIdx;"
    )

    # split into prompt and non-prompt
    # first do the leadingMuon
    if not skipLeading:
        rdfLeadingPrompt = rdf.Filter("leadingGenMuonIdx>-1").Filter("PromptMuonCond[leadingGenMuonIdx]")
        rdfLeadingNon = rdf.Filter("(leadingGenMuonIdx==-1) || (leadingGenMuonIdx>-1 && PromptMuonCond[leadingGenMuonIdx]==0)")
        rdfLeadingPrompt.Snapshot("Events", outLeadingPromptFile)
        rdfLeadingNon.Snapshot("Events", outLeadingNonFile)

    if not skipSubleading:
        rdfSubLeadingPrompt = rdf.Filter("subleadingGenMuonIdx>-1").Filter("PromptMuonCond[subleadingGenMuonIdx]")
        rdfSubLeadingNon = rdf.Filter("(subleadingGenMuonIdx==-1) || (subleadingGenMuonIdx>-1 && PromptMuonCond[subleadingGenMuonIdx]==0)")
        rdfSubLeadingPrompt.Snapshot("Events", outSubLeadingPromptFile)
        rdfSubLeadingNon.Snapshot("Events", outSubLeadingNonFile)
  

def main():
    this_dir = os.path.dirname(os.path.abspath(__file__))

    # specify the periods to run
    periods = [
        "RunIII2024Summer24NanoAODv15"
    ]

    # only developed for Run2024
    targetDS = [
        "ttbarDL", "ttbarSL", "ttbar4Q", 
        "TTBB_DL", "TTBB_SL", "TTBB_4Q",
        "TTW", "TTZ_low", "TTZ_high",
        "TTHBB", "TTHnonBB",
        "TTTT"
    ]

    dirPt1 = "/scratch/tiyang/S1_directly_SS/mc/RunIII2024Summer24NanoAODv15/"
    dsPt1 = []
    for ds in targetDS:
        if ds in os.listdir(dirPt1):
            dsPt1.append(ds)

    dirPt2 = "/scratch/tiyang/using2023For2024_dimuonSS/mc/RunIII2024Summer24NanoAODv15/"
    dsPt2 = []
    for ds in targetDS:
        if ds in os.listdir(dirPt2):
            dsPt2.append(ds)

    for era in periods:
        print("split era: " + str(era))
        for ds in dsPt1:
            filePaths = [dirPt1 + ds + "/" + i for i in os.listdir(dirPt1+ds)]
            split(era, ds, filePaths)
        for ds in dsPt2:
            filePaths = [dirPt2 + ds + "/" + i for i in os.listdir(dirPt2+ds)]
            split(era, ds, filePaths)

if __name__ == "__main__":
    main()
