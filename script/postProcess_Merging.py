import ROOT
import correctionlib
import json
import os
from abc import ABC, abstractmethod

class correctionApplier(ABC):
    # must keep the chain structure
    # different from c++, chain will be GCed and rdf implicitly rely on it
    def __init__(self, jsonPath, era, dataset):
        with open(jsonPath) as jFile:
            jsonFull = json.load(jFile)
        self.era = era
        self.dataset = dataset
        fileDir = jsonFull["dir"][dataset]
        fileList = jsonFull["file"][dataset]
        self.ifExist = 1
        self.c1 = ROOT.TChain("Events")
        self.numberOfFiles = 0 
        for fUnique in fileList:
            filePath = fileDir+"/"+fUnique
            if os.path.isfile(filePath):
                self.c1.Add(filePath)
                self.numberOfFiles += 1
        print("chain loaded")

    @abstractmethod
    def reweight_and_save(self):
        pass

def initializing():
    ROOT.EnableImplicitMT()
    correctionlib.register_pyroot_binding()
    ROOT.gInterpreter.Declare("""
        #include "correction.h"
        #include <memory>
        std::unique_ptr<correction::CorrectionSet> cset;

        void load_correction(const std::string& filepath) {
            cset = correction::CorrectionSet::from_file(filepath);
        }
    """)

def correction(sampleList, era, datasets, operationClass):
    if not issubclass(operationClass, correctionApplier):
        raise TypeError("operationClass must inherit from correctionApplier!")
    for dataset in datasets:
        print("skim dataset: " + str(dataset))
        eraProcessor = operationClass(sampleList[era], era, dataset)
        if eraProcessor.numberOfFiles == 0:
            continue
        eraProcessor.reweight_and_save()