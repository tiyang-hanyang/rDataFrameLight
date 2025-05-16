import postProcess
import ROOT
import correctionlib
import os

# redefine the branch operation method and directory to save
class ZPtProcesser(postProcess.correctionApplier):
    def reweight_and_save(self):
        # preparing the computation of JetVetoMap
        rdf = ROOT.ROOT.RDataFrame(self.c1)
        print("dataframe created")
        branchArray = list(rdf.GetColumnNames())
        branchArray.append("ZptWgt")

        rdf = rdf.Define("nGenZ", "GenPart_pdgId[(GenPart_pdgId == 23) && (GenPart_status==62)].size()")
        rdf = rdf.Define("Zpt", "nGenZ ? GenPart_pt[(GenPart_pdgId == 23) && (GenPart_status==62)][0] : 0.0")
        rdf = rdf.Define("ZptWgt",'nGenZ ? corr->evaluate({"LO", Zpt, "nom"}) : 1.0')
        print("ZptWgt DONE")

        # saving with corresponding names
        print("preparing save options")
        opt = ROOT.RDF.RSnapshotOptions()
        opt.fCompressionLevel = 9

        print("create out folder")
        outDir = "./ZPtCorrection/"+self.era+"/"
        if not os.path.exists(outDir):
            os.makedirs(outDir)
        outputPath = outDir+self.dataset+"_skimmed_ZPt.root"
        print("save")
        rdf.Snapshot("Events", outputPath, branchArray, opt)
        print("done")

def main():
    postProcess.initializing()
    this_dir = os.path.dirname(os.path.abspath(__file__))
    ROOT.gInterpreter.AddIncludePath(this_dir)
    ROOT.gInterpreter.ProcessLine('#include "ZPtcorrection.C"')

    # specify the periods to run
    periods = [
        "Run3Summer23NanoAODv12", 
        "Run3Summer23BPixNanoAODv12"
    ]

    # correctionlib json files accordingly
    correctionFiles = {
        "Run3Summer23NanoAODv12": "/home/tiyang/public/rDataFrameLight_git/correction/Zpt/DY_pTll_weights_2023preBPix_v2.json",
        "Run3Summer23BPixNanoAODv12": "/home/tiyang/public/rDataFrameLight_git/correction/Zpt/DY_pTll_weights_2023postBPix_v2.json"
    }

    from MetricSkimmedFiles import ZPtFileLists as jsonLists
    # from MetricSkimmedFiles import datasets
    datasets = {
        "Run3Summer23NanoAODv12": ["DY2L"],
        "Run3Summer23BPixNanoAODv12": ["DY2L"]
    }

    for era in periods:
        print("skim era: " + str(era))
        ROOT.gInterpreter.ProcessLine('load_correction("' + correctionFiles[era] +'");')
        ROOT.gInterpreter.ProcessLine('bind_correction();')
        postProcess.correction(jsonLists, era, datasets[era], ZPtProcesser)

if __name__ == "__main__":
    main()