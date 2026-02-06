import os
import sys
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, ".."))
from jobDef import GeneralJob


class RunIII2024Summer24NanoAODv15_job(GeneralJob):
    def __init__(self):
        super().__init__()
        self.periods = ["Run2024C"]
        self.datasets = [
            "Muon0",
            "Muon1"
        ]
        self.mergeDS = self.datasets

        # data corrections
        self.workflow = [
            # ("Rcorr", "RochesterCorr_MC.py"), 
            ("JEC", "JEC_Data.py"),
            ("JVM", "JetVetoMap.py"),
        ]

        self.fileJson = {
            "Run2024C": "../../../json/samples/special_relax_miniisolation/Run2024C_special_relax_miniIso.json",
        }

        for (era, fPath) in self.fileJson.items():
            self.fileJson[era] = os.path.abspath(os.path.join(here, fPath))


        self.outDir = "/home/tiyang/public/AR_special_relax_miniIsolation_0128_corr/"

    def declare(self):
        return self.periods, self.datasets, self.mergeDS, self.workflow, self.fileJson, self.outDir
