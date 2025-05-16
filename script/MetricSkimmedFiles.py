RawSampleLists = {
    # mc 
    "Run3Summer23NanoAODv12": "json/samples/Run3Summer23NanoAODv12.json",
    "Run3Summer23BPixNanoAODv12": "json/samples/Run3Summer23BPixNanoAODv12.json",

    # data
    "Run2023C": "json/samples/Run2023C.json",
    "Run2023D": "json/samples/Run2023D.json"
}

MetricSkimmedLists = {
    # mc
    "Run3Summer22EENanoAODv12": "json/samples/MetricMuon/MetricMuonSkimmed_Run3Summer22EENanoAODv12.json",
    "Run3Summer23NanoAODv12": "json/samples/MetricMuon/MetricMuonSkimmed_Run3Summer23NanoAODv12.json",
    "Run3Summer23BPixNanoAODv12": "json/samples/MetricMuon/MetricMuonSkimmed_Run3Summer23BPixNanoAODv12.json",

    # data
    "Run2023C": "json/samples/MetricMuon/MetricMuonSkimmed_Run2023C.json",
    "Run2023D": "json/samples/MetricMuon/MetricMuonSkimmed_Run2023D.json"
}

JetVetoedLists = {
    "Run3Summer23NanoAODv12": "json/samples/JetVeto/JetVeto_Run3Summer23NanoAODv12.json",
    "Run3Summer23BPixNanoAODv12": "json/samples/JetVeto/JetVeto_Run3Summer23BPixNanoAODv12.json",

    "Run2023C": "json/samples/JetVeto/JetVeto_Run2023C.json",
    "Run2023D": "json/samples/JetVeto/JetVeto_Run2023D.json"
}

RochesterCorrectedLists = {
    "Run3Summer23NanoAODv12": "json/samples/Rochester/Rochester_Run3Summer23NanoAODv12.json",
    "Run3Summer23BPixNanoAODv12": "json/samples/Rochester/Rochester_Run3Summer23BPixNanoAODv12.json"
}

PUFileLists = {
    "Run3Summer23NanoAODv12": "json/samples/Run3Summer23NanoAODv12_ForPUWeight.json", 
    "Run3Summer23BPixNanoAODv12": "json/samples/Run3Summer23BPixNanoAODv12_ForPUWeight.json"
}

MuonCorrectionFileLists = {
    "Run3Summer23NanoAODv12": "json/samples/Run3Summer23NanoAODv12_MuonCorr.json", 
    "Run3Summer23BPixNanoAODv12": "json/samples/Run3Summer23BPixNanoAODv12_MuonCorr.json"
}

MuonDS = ["Muon0", "Muon1"]
MuonDS22 = ["Muon"]

mcDatasets = [
    # DY like
    "DY2L",
    "WJets",

    # ttbar
    "ttbarDL",
    "ttbarSL",

    # VV
    "WW2L2Nu",
    "WZ2L2Q",
    "ZZ2L2Nu",
    "ZZ2L2Q",
    "WZ3l",
    "ZZ4l",

    # VVV
    "WWW",
    "WWZ",
    "WZZ",
    "ZZZ",

    # tV
    "TWm2L",
    "TbarWp2L",
    "TZQB",

    # ttV
    "TTZ_low",
    "TTZ_high",
    "TTW",

    # ttVV
    "TTWW",
    "TTWZ",
    "TTZZ",

    # 4t
    "TTTT",

    # Higgs
    "TTHBB",
    "TTHnonBB",
    "THQ",
    "THW",
    "TTWH",
    "TTZH",

    # QCD
    "QCD_15_20_mu",
    "QCD_20_30_mu",
    "QCD_30_50_mu",
    "QCD_50_80_mu",
    "QCD_80_120_mu",
    "QCD_120_170_mu",
    "QCD_170_300_mu",
    "QCD_300_470_mu",
    "QCD_470_600_mu",
    "QCD_600_800_mu",
    "QCD_800_1000_mu",
    "QCD_1000_mu",

    # signal
    "TTHH_DL_2B2V",
    "TTHH_SL_2B2V",

    # gluons
    "GGZZ_2mu2e",
    "GGZZ_4mu",
    "GGZZ_4e",
    "GGZZ_2mu2tau",
    "GGZZ_2e2tau",
    "DYG",
    "WG"
]

# to import from this
datasets = {
    "2022": MuonDS22,
    "Run2023C": MuonDS,
    "Run2023D": MuonDS,
    "Run2024C": MuonDS,
    "Run3Summer22EENanoAODv12": mcDatasets,
    "Run3Summer23NanoAODv12": mcDatasets,
    "Run3Summer23BPixNanoAODv12": mcDatasets,
}