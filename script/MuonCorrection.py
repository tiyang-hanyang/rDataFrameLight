import postProcess
import ROOT
import correctionlib
import os

# redefine the branch operation method and directory to save
class MuonScaleProcesser(postProcess.correctionApplier):
    def reweight_and_save(self):
        # preparing the computation of JetVetoMap
        rdf = ROOT.ROOT.RDataFrame(self.c1)
        print("dataframe created")
        # branchArray = list(rdf.GetColumnNames())
        # branchArray.append("Muon_IDscale")
        # branchArray.append("Muon_Isoscale")
        branchArray = ROOT.std.vector("string")()
        for name in rdf.GetColumnNames():
            branchArray.push_back(name)
        branchArray.push_back("Muon_IDscale")
        branchArray.push_back("Muon_Isoscale")
        branchArray.push_back("MuonScale")

        rdf = rdf.Filter("(MuonPresel == 1)")
        rdf = rdf.Define("Muon_IDscale", "MuonIDScale(Muon_eta, Muon_pt_Rcorr)")
        rdf = rdf.Define("Muon_Isoscale", "MuonIsoScale(Muon_eta, Muon_pt_Rcorr)")
        rdf = rdf.Define("MuonScale", "EventMuonSF(Muon_IDscale, Muon_Isoscale, isGoodMuon_25, isGoodMuon_15)")

        print("Muon correction DONE")

        # saving with corresponding names
        print("preparing save options")
        opt = ROOT.RDF.RSnapshotOptions()
        opt.fCompressionLevel = 9

        print("create out folder")
        outDir = "./DiMuonCorrected/"+self.era+"/"
        if not os.path.exists(outDir):
            os.makedirs(outDir)
        outputPath = outDir+self.dataset+"_skimmed_DiMuonCorrected.root"
        print("save")
        rdf.Snapshot("Events", outputPath, branchArray, opt)
        print("done")

def main():
    postProcess.initializing()
    this_dir = os.path.dirname(os.path.abspath(__file__))
    ROOT.gInterpreter.AddIncludePath(this_dir)
    ROOT.gInterpreter.ProcessLine('#include "MuonCorrection.C"')

    # specify the periods to run
    periods = [
        "Run3Summer23NanoAODv12", 
        # Run2023D not ready yet
        # "Run3Summer23BPixNanoAODv12"
    ]

    # correctionlib json files accordingly
    correctionFiles = {
        # "Run3Summer22EENanoAODv12": "/home/tiyang/public/rDataFrameLight_git/correction/",
        "Run3Summer23NanoAODv12": this_dir+ "/../../" +"correction/POGCorr/POG/MUO/2023_Summer23/muon_Z.json",
        "Run3Summer23BPixNanoAODv12": this_dir+ "/../../" +"correction/POGCorr/POG/MUO/2023_Summer23/muon_Z.json"
    }

    from MetricSkimmedFiles import MuonCorrectionFileLists as jsonLists
        for key in jsonLists.keys():
        jsonLists[key] = this_dir+ "/../" + jsonLists[key]

    # from MetricSkimmedFiles import datasets
    datasets = ["DY2L"]

    for era in periods:
        print("skim era: " + str(era))
        ROOT.gInterpreter.ProcessLine('load_correction("' + correctionFiles[era] +'");')
        # ROOT.gInterpreter.ProcessLine('bind_correction("' + era + '");')
        postProcess.correction(jsonLists, era, datasets, MuonScaleProcesser)

if __name__ == "__main__":
    main()