import postProcess
import ROOT
import correctionlib
import os

# redefine the branch operation method and directory to save
class PUReweightProcesser(postProcess.correctionApplier):
    def reweight_and_save(self):
        # preparing the computation of JetVetoMap
        rdf = ROOT.ROOT.RDataFrame(self.c1)
        print("dataframe created")
        branchArray = list(rdf.GetColumnNames())
        branchArray.append("PUWeight")

#############
        rdf = rdf.Define("PUWeight", "PUReweightFunc(Pileup_nTrueInt)")

        print("PU reweighting DONE")

        # saving with corresponding names
        print("preparing save options")
        opt = ROOT.RDF.RSnapshotOptions()
        opt.fCompressionLevel = 9

        print("create out folder")
        outDir = "./PUWeight/"+self.era+"/"
        if not os.path.exists(outDir):
            os.makedirs(outDir)
        outputPath = outDir+self.dataset+"_skimmed_PUWeight.root"
        print("save")
        rdf.Snapshot("Events", outputPath, branchArray, opt)
        print("done")

def main():
    postProcess.initializing()
    this_dir = os.path.dirname(os.path.abspath(__file__))
    ROOT.gInterpreter.AddIncludePath(this_dir)
    ROOT.gInterpreter.ProcessLine('#include "PUWeight.C"')

    # specify the periods to run
    periods = [
        "Run3Summer23NanoAODv12", 
        # Run2023D not ready yet
        # "Run3Summer23BPixNanoAODv12"
    ]

    # correctionlib json files accordingly
    correctionFiles = {
        # "Run3Summer22EENanoAODv12": "/home/tiyang/public/rDataFrameLight_git/correction/",
        "Run3Summer23NanoAODv12": "/home/tiyang/public/rDataFrameLight_git/correction/jsonpog-integration-master/POG/LUM/2023_Summer23/puWeights.json",
        "Run3Summer23BPixNanoAODv12": "/home/tiyang/public/rDataFrameLight_git/correction/jsonpog-integration-master/POG/LUM/2023_Summer23BPix/puWeights.json"
    }

    from MetricSkimmedFiles import PUFileLists as jsonLists
    # from MetricSkimmedFiles import datasets

    # firstly only test the DY2L
    datasets = ["DY2L"]

    for era in periods:
        print("skim era: " + str(era))
        ROOT.gInterpreter.ProcessLine('load_correction("' + correctionFiles[era] +'");')
        ROOT.gInterpreter.ProcessLine('bind_correction("' + era + '");')
        # postProcess.correction(jsonLists, era, datasets[era], PUReweightProcesser)
        postProcess.correction(jsonLists, era, datasets, PUReweightProcesser)

if __name__ == "__main__":
    main()