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
            ("JVM", "JetVetoMap_AR_rough.py"),
        ]

        self.fileJson = {
            "Run2024C": "../../../json/samples/applicationRegion/Run2024C_roughSel.json",
            "Run2024D": "../../../json/samples/applicationRegion/Run2024D_roughSel.json",
            "Run2024E": "../../../json/samples/applicationRegion/Run2024E_roughSel.json",
            "Run2024F": "../../../json/samples/applicationRegion/Run2024F_roughSel.json",
            "Run2024G": "../../../json/samples/applicationRegion/Run2024G_roughSel.json",
            "Run2024H": "../../../json/samples/applicationRegion/Run2024H_roughSel.json",
            "Run2024I": "../../../json/samples/applicationRegion/Run2024I_roughSel.json",
        }

        for (era, fPath) in self.fileJson.items():
            self.fileJson[era] = os.path.abspath(os.path.join(here, fPath))


        self.outDir = "/home/tiyang/public/AR_0204_corr/"

    def declare(self):
        return self.periods, self.datasets, self.mergeDS, self.workflow, self.fileJson, self.outDir
