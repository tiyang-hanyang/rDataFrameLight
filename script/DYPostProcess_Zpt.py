import ROOT
import correctionlib
import json
import os

class DYZptWeighter:
    # must keep the chain structure
    # different from c++, chain will be GCed and rdf implicitly rely on it
    def __init__(self, jsonPath, era):
        with open(jsonPath) as jFile:
            jsonFull = json.load(jFile)
        self.era = era

        DYDir = jsonFull["dir"]["DY2L"]
        DYList = jsonFull["file"]["DY2L"]
        self.c1 = ROOT.TChain("Events")
        for fPath in DYList:
            self.c1.Add(DYDir+"/"+fPath)
        print("chain loaded")

    def reweight_and_save(self):
        rdf = ROOT.ROOT.RDataFrame(self.c1)
        print("dataframe created")
        # only record the ZptWgt in the definition
        branchArray = list(rdf.GetColumnNames())
        branchArray.append("ZptWgt")

        # do the definition
        rdf = rdf.Define("nGenZ", "GenPart_pdgId[(GenPart_pdgId == 23) && (GenPart_status==62)].size()")
        rdf = rdf.Define("Zpt", "nGenZ ? GenPart_pt[(GenPart_pdgId == 23) && (GenPart_status==62)][0] : 0.0")
        rdf = rdf.Define("ZptWgt",'nGenZ ? corr->evaluate({"LO", Zpt, "nom"}) : 1.0')
        print("ZptWgt DONE")

        print("preparing save options")
        opt = ROOT.RDF.RSnapshotOptions()
        opt.fCompressionLevel = 9

        print("create out folder")
        outDir = "./DY2L_ZptWgt_added/"+self.era+"/"
        if not os.path.exists(outDir):
            os.makedirs(outDir)
        outputPath = outDir+"DY2L_skimmed_ZptWgt.root"
        print("save")
        rdf.Snapshot("Events", outputPath, branchArray, opt)
        print("done")

def main():
    ROOT.EnableImplicitMT()
    periods = ["Run2023C", "Run2023D"]
    jsonLists = {
        "Run2023C": "/home/tiyang/public/rDataFrameLight_git/source/json/samples/MetricMuon/MetricMuonSkimmed_Run2023C.json",
        "Run2023D": "/home/tiyang/public/rDataFrameLight_git/source/json/samples/MetricMuon/MetricMuonSkimmed_Run2023D.json"
    }
    correctionlib.register_pyroot_binding()
    ROOT.gInterpreter.Declare('auto cset = correction::CorrectionSet::from_file("/home/tiyang/public/rDataFrameLight_git/correction/DY_pTll_weights_2023postBPix_v2.json");')
    ROOT.gInterpreter.Declare('auto corr = cset->at("DY_pTll_reweighting");')
    for era in periods:
        eraProcessor = DYZptWeighter(jsonLists[era], era)
        eraProcessor.reweight_and_save()

if __name__ == "__main__":
    main()