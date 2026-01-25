import ROOT

def checkPromptInfo(era, dataset):
    files = [
        "splitPrompt/"+era+"/"+dataset+"/"+dataset+"_leadingMuonPrompt_incluTau.root",
        "splitPrompt/"+era+"/"+dataset+"/"+dataset+"_leadingMuonNonPrompt_incluTau.root",
    ]

    rdf = ROOT.RDataFrame("Events", files[0])
    rdf = rdf.Define("nTruthPrompt", "Nonzero(PromptMuonCond).size()")
    nLeadingPrompt = rdf.Sum("weight_XS").GetValue() * 108.96
    nZeroGenPrompt = rdf.Filter("nTruthPrompt==0").Sum("weight_XS").GetValue() * 108.96
    nOneGenPrompt = rdf.Filter("nTruthPrompt==1").Sum("weight_XS").GetValue() * 108.96
    nTwoGenPrompt = rdf.Filter("nTruthPrompt==2").Sum("weight_XS").GetValue() * 108.96
    nMoreGenPrompt = rdf.Filter("nTruthPrompt>2").Sum("weight_XS").GetValue() * 108.96
    nTwoRecoPrompt = rdf.Filter("subleadingGenMuonIdx>-1").Filter("PromptMuonCond[subleadingGenMuonIdx]").Sum("weight_XS").GetValue() * 108.96
    nOneRecoPrompt = nLeadingPrompt-nTwoRecoPrompt

    rdf = ROOT.RDataFrame("Events", files[1])
    rdf = rdf.Define("nTruthPrompt", "Nonzero(PromptMuonCond).size()")
    nLeadingNonPrompt = rdf.Sum("weight_XS").GetValue() * 108.96
    nZeroGenPrompt += rdf.Filter("nTruthPrompt==0").Sum("weight_XS").GetValue() * 108.96
    nOneGenPrompt += rdf.Filter("nTruthPrompt==1").Sum("weight_XS").GetValue() * 108.96
    nTwoGenPrompt += rdf.Filter("nTruthPrompt==2").Sum("weight_XS").GetValue() * 108.96
    nMoreGenPrompt += rdf.Filter("nTruthPrompt>2").Sum("weight_XS").GetValue() * 108.96
    nSubleadingPrompt = rdf.Filter("subleadingGenMuonIdx>-1").Filter("PromptMuonCond[subleadingGenMuonIdx]").Sum("weight_XS").GetValue() * 108.96
    nZeroRecoPrompt = nLeadingNonPrompt - nSubleadingPrompt
    nOneRecoPrompt += nSubleadingPrompt

    print(dataset, "reco prompt count: 0, 1, 2")
    print("\t", nZeroRecoPrompt, ",", nOneRecoPrompt, ",", nTwoRecoPrompt)
    print(dataset, "genpart prompt count: 0, 1, 2, 3-")
    print("\t", nZeroGenPrompt, ",", nOneGenPrompt, ",", nTwoGenPrompt, ",", nMoreGenPrompt)


def main():
    # specify the periods to run
    periods = [
        "RunIII2024Summer24NanoAODv15"
    ]

    # only developed for Run2024
    targetDS = [
        "ttbarDL", "ttbarSL",
        "TTBB_DL", "TTBB_SL",
        "TTW", "TTZ_low", "TTZ_high",
        "TTHBB", "TTHnonBB",
        "TTTT"
    ]

    for era in periods:
        for ds in targetDS:
            checkPromptInfo(era, ds)

if __name__ == "__main__":
    main()
