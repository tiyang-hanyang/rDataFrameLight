import os
import sys
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, ".."))
from jobDef import GeneralJob

class RunIII2024Summer24NanoAODv15_job(GeneralJob):
    def __init__(self):
        super().__init__()
        self.periods = ["RunIII2024Summer24NanoAODv15_AR"]
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

            # VV
            "WW",
            "WZ",
            "ZZ",

            # WJet
            "WJet_1J",
            "WJet_2J",
            "WJet_3J",
            "WJet_4J"
        ]
        self.mergeDS = self.datasets

        self.workflow = [
            # b-tag correction after the jet selections
            ("BTag", "BTagCorr.py"),
        ]

        self.fileJson = {
            "RunIII2024Summer24NanoAODv15_AR": "../../../json/samples/applicationRegion/RunIII2024Summer24NanoAODv15_fullFakeableSel.json",
        }
        for (era, fPath) in self.fileJson.items():
            self.fileJson[era] = os.path.abspath(os.path.join(here, fPath))

        self.outDir = "/home/tiyang/public/AR_0204_btag_corrected/"

    def declare(self):
        return self.periods, self.datasets, self.mergeDS, self.workflow, self.fileJson, self.outDir
