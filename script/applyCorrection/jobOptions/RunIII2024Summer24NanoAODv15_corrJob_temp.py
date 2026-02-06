import os
import sys
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, ".."))
from jobDef import GeneralJob

class RunIII2024Summer24NanoAODv15_job(GeneralJob):
    def __init__(self):
        super().__init__()
        self.periods = ["RunIII2024Summer24NanoAODv15"]
        self.datasets = [
            # DY
            "DY2Mu",
            "DY2Mu_low",

            # ttbar
            "ttbarDL",
            "ttbarSL",
            "ttbar4Q",

            # ttbb
            "TTBB_DL",
            "TTBB_SL",
            "TTBB_4Q",

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

            # ttH
            "TTHBB",
            "TTHnonBB",

            # tW
            "TWm2L",
            "TbarWp2L",
            "TWm1L",
            "TbarWp1L",
            "TZQB",

            # VV
            "WW2L2Nu",
            "WZ2L2Q",
            "ZZ2L2Q",
            "ZZ2L2Nu",
            "WZ3l",
            "ZZ4l",

            "WW",
            "WZ",
            "ZZ",

            # VVV
            "WWW",
            "WWZ",
            "WZZ",
            "ZZZ",

            # WJet
            "WJet_1J",
            "WJet_2J",
            "WJet_3J",
            "WJet_4J",

            # TTZ
            "TTZ_low",
            "TTZ_high",

            # signal
            # "TTHH_DL_2B2W_batch1",
            # "TTHH_SL_2B2W_batch1",
        ]
        self.mergeDS = self.datasets

        # here is one subtle thing, the btag needs the ds-dependent tagging efficiency instead of the only weight
        # this makes it cannot only extracted the weights directly
        self.workflow = [
            # ("Rcorr", "RochesterCorr_MC.py"), 
            ("JEC", "JEC_MC.py"),
            ("JVM", "JetVetoMap.py"),
            ("PU", "PUWeight.py"), 
            # muon correction need to also define what is good muon, and only correct the good muons
            # ("Muon", "Muon_corr.py"),
            # the b-tagging validation need reliable Jet correction, thus not do yet
            # ("BTag", "BTagCorr.py"),
        ]

        self.fileJson = {
            "RunIII2024Summer24NanoAODv15": "../../../json/samples/Dimuon_NanoAOD/RunIII2024Summer24NanoAODv15_forCorr_temp.json",
        }

        for (era, fPath) in self.fileJson.items():
            self.fileJson[era] = os.path.abspath(os.path.join(here, fPath))

        self.outDir = "/home/tiyang/public/SR_medium_muon_0127/"

    def declare(self):
        return self.periods, self.datasets, self.mergeDS, self.workflow, self.fileJson, self.outDir
