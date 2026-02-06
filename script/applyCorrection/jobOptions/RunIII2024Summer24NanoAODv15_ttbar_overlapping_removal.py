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
            # ttbar
            "ttbarDL",
            "ttbarSL",
            "ttbar4Q"
        ]
        self.mergeDS = self.datasets

        # here is one subtle thing, the btag needs the ds-dependent tagging efficiency instead of the only weight
        # this makes it cannot only extracted the weights directly
        self.workflow = [
            ("overlappingRemoval", "ttbar_overlapping_removal.py"),
        ]

        # one possible thing. Maybe I should not have the overlapping check after the dimuon matching, but instead should check at the very origin level? 
        # whatever, finally this ttbar will be used with the dimuon cut, but the ratio might be more clear without the selection.
        self.fileJson = {
            "RunIII2024Summer24NanoAODv15": "../../../json/samples/Dimuon_NanoAOD/PUJECcorrected/RunIII2024Summer24NanoAODv15_corrected_temp_pre_ttbarRemoval.json",
        }
        for (era, fPath) in self.fileJson.items():
            self.fileJson[era] = os.path.abspath(os.path.join(here, fPath))

        self.outDir = "/home/tiyang/public/SR_medium_muon_0127_ttbar_overlapping_removal/"

    def declare(self):
        return self.periods, self.datasets, self.mergeDS, self.workflow, self.fileJson, self.outDir