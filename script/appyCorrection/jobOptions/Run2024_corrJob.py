import os
import sys
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, ".."))
from jobDef import GeneralJob


class RunIII2024Summer24NanoAODv15_job(GeneralJob):
    def __init__(self):
        super().__init__()
        self.periods = ["Run2024C", "Run2024D", "Run2024E", "Run2024F", "Run2024G", "Run2024H", "Run2024I"]
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
            "Run2024C": "../../../json/samples/Dimuon_NanoAOD/Run2024C_forCorr.json",
            "Run2024D": "../../../json/samples/Dimuon_NanoAOD/Run2024D_forCorr.json",
            "Run2024E": "../../../json/samples/Dimuon_NanoAOD/Run2024E_forCorr.json",
            "Run2024F": "../../../json/samples/Dimuon_NanoAOD/Run2024F_forCorr.json",
            "Run2024G": "../../../json/samples/Dimuon_NanoAOD/Run2024G_forCorr.json",
            "Run2024H": "../../../json/samples/Dimuon_NanoAOD/Run2024H_forCorr.json",
            "Run2024I": "../../../json/samples/Dimuon_NanoAOD/Run2024I_forCorr.json",
        }

        for (era, fPath) in self.fileJson.items():
            self.fileJson[era] = os.path.abspath(os.path.join(here, fPath))


        self.outDir = "/home/tiyang/public/SR_medium_muon_0127/"

    def declare(self):
        return self.periods, self.datasets, self.mergeDS, self.workflow, self.fileJson, self.outDir