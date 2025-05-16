import postProcess
import ROOT
import correctionlib
import os

# redefine the branch operation method and directory to save
class JVMProcesser(postProcess.correctionApplier):
    def reweight_and_save(self):
        # preparing the computation of JetVetoMap
        rdf = ROOT.ROOT.RDataFrame(self.c1)
        print("dataframe created")
        branchArray = list(rdf.GetColumnNames())
        branchArray.append("PassVetoJet_jetId")
        branchArray.append("PassVetoJet_pt")
        branchArray.append("PassVetoJet_phi")
        branchArray.append("PassVetoJet_eta")
        branchArray.append("PassVetoJet_mass")
        branchArray.append("PassVetoJet_area")
        branchArray.append("nPassVetoJet")
        branchArray.append("PassVetoBJet_jetId")
        branchArray.append("PassVetoBJet_pt")
        branchArray.append("PassVetoBJet_phi")
        branchArray.append("PassVetoBJet_eta")
        branchArray.append("PassVetoBJet_mass")
        branchArray.append("PassVetoBJet_area")
        branchArray.append("nPassVetoBJet")

        rdf = rdf.Define("PassJetVeto", "passJetVetoFunc(GoodJet_eta, GoodJet_phi)")
        rdf = rdf.Define("PassVetoJet_jetId", "GoodJet_jetId[PassJetVeto]")
        rdf = rdf.Define("PassVetoJet_pt", "GoodJet_pt[PassJetVeto]")
        rdf = rdf.Define("PassVetoJet_eta", "GoodJet_eta[PassJetVeto]")
        rdf = rdf.Define("PassVetoJet_phi", "GoodJet_phi[PassJetVeto]")
        rdf = rdf.Define("PassVetoJet_mass", "GoodJet_mass[PassJetVeto]")
        rdf = rdf.Define("PassVetoJet_area", "GoodJet_area[PassJetVeto]")
        rdf = rdf.Define("nPassVetoJet", "PassVetoJet_jetId.size()")

        rdf = rdf.Define("PassJetVetoBJet", "passJetVetoFunc(BJet_eta, BJet_phi)")
        rdf = rdf.Define("PassVetoBJet_jetId", "BJet_jetId[PassJetVetoBJet]")
        rdf = rdf.Define("PassVetoBJet_pt", "BJet_pt[PassJetVetoBJet]")
        rdf = rdf.Define("PassVetoBJet_eta", "BJet_eta[PassJetVetoBJet]")
        rdf = rdf.Define("PassVetoBJet_phi", "BJet_phi[PassJetVetoBJet]")
        rdf = rdf.Define("PassVetoBJet_mass", "BJet_mass[PassJetVetoBJet]")
        rdf = rdf.Define("PassVetoBJet_area", "BJet_area[PassJetVetoBJet]")
        rdf = rdf.Define("nPassVetoBJet", "PassVetoBJet_jetId.size()")
        print("jetveomap DONE")

        # saving with corresponding names
        print("preparing save options")
        opt = ROOT.RDF.RSnapshotOptions()
        opt.fCompressionLevel = 9

        print("create out folder")
        outDir = "./passJetVeto/"+self.era+"/"
        if not os.path.exists(outDir):
            os.makedirs(outDir)
        outputPath = outDir+self.dataset+"_skimmed_jetvetomap.root"
        print("save")
        rdf.Snapshot("Events", outputPath, branchArray, opt)
        print("done")

def main():
    postProcess.initializing()
    this_dir = os.path.dirname(os.path.abspath(__file__))
    ROOT.gInterpreter.AddIncludePath(this_dir)
    ROOT.gInterpreter.ProcessLine('#include "JetVetoMap.C"')

    # specify the periods to run
    periods = [
        #"Run2023C",
        #"Run2023D",
        #"Run3Summer22EENanoAODv12",
        "Run3Summer23NanoAODv12", 
        "Run3Summer23BPixNanoAODv12"
    ]

    # correctionlib json files accordingly
    correctionFiles = {
        "Run3Summer22EENanoAODv12": "/home/tiyang/public/rDataFrameLight_git/correction/JEC/jsonpog-integration-master/POG/JME/2022_Summer22EE/jetvetomaps.json",
        "Run2023C": "/home/tiyang/public/rDataFrameLight_git/correction/JEC/jsonpog-integration-master/POG/JME/2023_Summer23/jetvetomaps.json",
        "Run3Summer23NanoAODv12": "/home/tiyang/public/rDataFrameLight_git/correction/JEC/jsonpog-integration-master/POG/JME/2023_Summer23/jetvetomaps.json",
        "Run2023D": "/home/tiyang/public/rDataFrameLight_git/correction/JEC/jsonpog-integration-master/POG/JME/2023_Summer23BPix/jetvetomaps.json",
        "Run3Summer23BPixNanoAODv12": "/home/tiyang/public/rDataFrameLight_git/correction/JEC/jsonpog-integration-master/POG/JME/2023_Summer23BPix/jetvetomaps.json"
    }

    from MetricSkimmedFiles import JetVetoFileLists as jsonLists
    from MetricSkimmedFiles import datasets

    for era in periods:
        print("skim era: " + str(era))
        ROOT.gInterpreter.ProcessLine('load_correction("' + correctionFiles[era] +'");')
        ROOT.gInterpreter.ProcessLine('bind_correction("' + era + '");')
        postProcess.correction(jsonLists, era, datasets[era], JVMProcesser)

if __name__ == "__main__":
    main()
