import os
import sys
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, ".."))
from jobDef import GeneralJob

class RunIII2024Summer24NanoAODv15_job(GeneralJob):
    def __init__(self):
        super().__init__()
        self.periods = ["Run3Summer23BPixNanoAODv12"]
        self.datasets = [
            "TTW",
            "TTTT"
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
            "Run3Summer23BPixNanoAODv12": "../../../json/samples/Dimuon_NanoAOD/Run3Summer23BPixNanoAODv12_forCorr_temp.json",
        }

        for (era, fPath) in self.fileJson.items():
            self.fileJson[era] = os.path.abspath(os.path.join(here, fPath))

        self.outDir = "/home/tiyang/public/SR_medium_muon_0127/"

    def declare(self):
        return self.periods, self.datasets, self.mergeDS, self.workflow, self.fileJson, self.outDir
