import os
import sys
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, ".."))
from jobDef import GeneralJob

class Run3Summer23BPixNanoAODv12_job(GeneralJob):
    def __init__(self):
        super().__init__()
        self.periods = ["Run3Summer23BPixNanoAODv12_SSCR"]
        self.datasets = [
            "TTW",
            "TTTT"
        ]
        self.mergeDS = self.datasets

        self.workflow = [
            # b-tag correction after the jet selections
            ("BTag", "BTagCorr.py"),
        ]

        self.fileJson = {
            "Run3Summer23BPixNanoAODv12_SSCR": "../../../json/samples/SameSign_CR/Run3Summer23BPixNanoAODv12_SSCR.json",
        }
        for (era, fPath) in self.fileJson.items():
            self.fileJson[era] = os.path.abspath(os.path.join(here, fPath))

        self.outDir = "/home/tiyang/public/SameSign_CR_0127_fourJet_btag_corrected/"

    def declare(self):
        return self.periods, self.datasets, self.mergeDS, self.workflow, self.fileJson, self.outDir
