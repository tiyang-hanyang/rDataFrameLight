import postProcess_Merging
import ROOT
import os

# redefine the branch operation method and directory to save
class JVMProcesser(postProcess_Merging.correctionApplier):
    def reweight_and_save(self):
        # preparing the computation of JetVetoMap
        rdf = ROOT.ROOT.RDataFrame(self.c1)
        print("dataframe created")
        branchArray = list(rdf.GetColumnNames())
        branchArray.append("PassJetVeto")

        rdf = rdf.Define("PassJetVeto", "passJetVetoFunc(Jet_eta, Jet_phi)")

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
    postProcess_Merging.initializing()
    this_dir = os.path.dirname(os.path.abspath(__file__))
    ROOT.gInterpreter.AddIncludePath(this_dir)
    ROOT.gInterpreter.ProcessLine('#include "JetVetoMap.C"')

    # specify the periods to run
    periods = [
        #"Run2022C",
        #"Run2022D",
        # "Run3Summer22NanoAODv12",
        #"Run2022E",
        #"Run2022F",
        #"Run2022G",
        # "Run3Summer22EENanoAODv12",
        #"Run2023C",
        # "Run3Summer23NanoAODv12", 
        #"Run2023D",
        # "Run3Summer23BPixNanoAODv12",
        #"Run2024C",
        #"Run2024D",
        #"Run2024E",
        #"Run2024F",
        #"Run2024G",
        #"Run2024H",
        #"Run2024I",
        "RunIII2024Summer24NanoAODv15"
    ]

    # correctionlib json files accordingly
    correctionFiles = {
        "Run2022C": "/home/tiyang/public/rDataFrameLight_git/correction/JEC/jsonpog-integration-master/POG/JME/2022_Summer22/jetvetomaps.json",
        "Run2022D": "/home/tiyang/public/rDataFrameLight_git/correction/JEC/jsonpog-integration-master/POG/JME/2022_Summer22/jetvetomaps.json",
        "Run3Summer22NanoAODv12": "/home/tiyang/public/rDataFrameLight_git/correction/JEC/jsonpog-integration-master/POG/JME/2022_Summer22/jetvetomaps.json",
        "Run2022E": "/home/tiyang/public/rDataFrameLight_git/correction/JEC/jsonpog-integration-master/POG/JME/2022_Summer22EE/jetvetomaps.json",
        "Run2022F": "/home/tiyang/public/rDataFrameLight_git/correction/JEC/jsonpog-integration-master/POG/JME/2022_Summer22EE/jetvetomaps.json",
        "Run2022G": "/home/tiyang/public/rDataFrameLight_git/correction/JEC/jsonpog-integration-master/POG/JME/2022_Summer22EE/jetvetomaps.json",
        "Run3Summer22EENanoAODv12": "/home/tiyang/public/rDataFrameLight_git/correction/JEC/jsonpog-integration-master/POG/JME/2022_Summer22EE/jetvetomaps.json",
        "Run2023C": "/home/tiyang/public/rDataFrameLight_git/correction/JEC/jsonpog-integration-master/POG/JME/2023_Summer23/jetvetomaps.json",
        "Run3Summer23NanoAODv12": "/home/tiyang/public/rDataFrameLight_git/correction/JEC/jsonpog-integration-master/POG/JME/2023_Summer23/jetvetomaps.json",
        "Run2023D": "/home/tiyang/public/rDataFrameLight_git/correction/JEC/jsonpog-integration-master/POG/JME/2023_Summer23BPix/jetvetomaps.json",
        "Run3Summer23BPixNanoAODv12": "/home/tiyang/public/rDataFrameLight_git/correction/JEC/jsonpog-integration-master/POG/JME/2023_Summer23BPix/jetvetomaps.json",
        "Run2024C": this_dir+ "/../../" +"correction/POGCorr/POG/JME/2024_Summer24/jetvetomaps.json",
        "Run2024D": this_dir+ "/../../" +"correction/POGCorr/POG/JME/2024_Summer24/jetvetomaps.json",
        "Run2024E": this_dir+ "/../../" +"correction/POGCorr/POG/JME/2024_Summer24/jetvetomaps.json",
        "Run2024F": this_dir+ "/../../" +"correction/POGCorr/POG/JME/2024_Summer24/jetvetomaps.json",
        "Run2024G": this_dir+ "/../../" +"correction/POGCorr/POG/JME/2024_Summer24/jetvetomaps.json",
        "Run2024H": this_dir+ "/../../" +"correction/POGCorr/POG/JME/2024_Summer24/jetvetomaps.json",
        "Run2024I": this_dir+ "/../../" +"correction/POGCorr/POG/JME/2024_Summer24/jetvetomaps.json",
        "RunIII2024Summer24NanoAODv15": this_dir+ "/../../" +"correction/POGCorr/POG/JME/2024_Summer24/jetvetomaps.json"
    }

    from MetricSkimmedFiles import JetVetoFileLists as jsonLists
    for key in jsonLists.keys():
        jsonLists[key] = this_dir+ "/../" + jsonLists[key]

    # from MetricSkimmedFiles import datasets
    datasets = [
        "DY2L",
        "DY2L_low",
        "ttbarDL",
        "ttbarSL",
        "ttbar4Q",
        "QCD_80_120_mu",
        "QCD_120_170_mu",
        "QCD_170_300_mu",
        "TbarWp2L",
        "TWm2L",
        "TTHBB",
        "TTHnonBB",
    ]
    # datasets = ["Muon0", "Muon1"]

    for era in periods:
        print("skim era: " + str(era))
        ROOT.gInterpreter.ProcessLine('load_correction("' + correctionFiles[era] +'");')
        ROOT.gInterpreter.ProcessLine('bind_correction("' + era + '");')
        # postProcess.correction(jsonLists, era, datasets[era], JVMProcesser)
        postProcess_Merging.correction(jsonLists, era, datasets, JVMProcesser)

if __name__ == "__main__":
    main()
