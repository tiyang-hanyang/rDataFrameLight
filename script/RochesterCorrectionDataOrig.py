import postProcess
import ROOT
import os
import json

# redefine the branch operation method and directory to save
class RochesterCorrection(postProcess.correctionApplier):
    def reweight_and_save(self):
        for filePath in self.filePaths:
            # preparing the computation of JetVetoMap
            rdf = ROOT.RDataFrame("Events", filePath)
            print("dataframe loaded for:", filePath)
            existBranchArray = list(rdf.GetColumnNames())

            # get original wanted branch list
            this_dir = os.path.dirname(os.path.abspath(__file__))
            with open(this_dir+ "/../" +"json/branches/MuonMetricBranchNew.json") as branchJson:
                branchArray = json.load(branchJson)
            branchArray = [brname for brname in branchArray if brname in existBranchArray]

            # Rochester correction added 
            # data case, only correct pt scale, name to be synchronized with the final correction from mc
            branchArray.append("Muon_pt_Rcorr")
            rdf = rdf.Define('Muon_pt_Rcorr', 'pt_scale(1, Muon_pt, Muon_eta, Muon_phi, Muon_charge)')

            # saving with corresponding names
            print("preparing save options")
            opt = ROOT.RDF.RSnapshotOptions()
            opt.fCompressionLevel = 9

            print("create out folder")
            outDir = "/data2/common/skimmed_NanoAOD/skim_1024_RochesterCorr/data/"+self.era+"/"+self.dataset+"/"
            if not os.path.exists(outDir):
                os.makedirs(outDir)
            
            newFilePath = "-".join(filePath.split("/")[-2:])[:-5]+"_Rcorr.root"
            outputPath = outDir + newFilePath
            print("save")
            rdf.Snapshot("Events", outputPath, branchArray, opt)
            print("done")

def main():
    postProcess.initializing()
    this_dir = os.path.dirname(os.path.abspath(__file__))
    ROOT.gInterpreter.AddIncludePath(this_dir)
    ROOT.gInterpreter.ProcessLine('#include "MuonScaReOrig.cc"')

    # specify the periods to run
    periods = [
        #"Run2023C", 
        #"Run2023D",
        #"Run2022C",
        #"Run2022D",
        #"Run2022E",
        #"Run2022F",
        #"Run2022G",
        #"Run2024C",
        #"Run2024D",
        #"Run2024E",
        #"Run2024F",
        #"Run2024G",
        "Run2024H",
        "Run2024I"
    ]

    # correctionlib json files accordingly
    # documentation: https://cms-analysis-corrections.docs.cern.ch/corrections/MUO/Run3-23DSep23-Summer23BPix-NanoAODv12/latest/#muon_zjsongz
    # 2022 and 2023 already merged into the POG documentation https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/MUO
    # previously reading from https://gitlab.cern.ch/cms-muonPOG/muonscarekit/-/tree/master?ref_type=heads
    # correction method also from this page
    # 2024 from https://gitlab.cern.ch/cms-analysis-corrections/MUO/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15, not in the POG json collection yet
    correctionFiles = {
        "Run2023C": this_dir+ "/../../" +"correction/POGCorr/POG/MUO/2023_Summer23/muon_scalesmearing.json",
        "Run2023D": this_dir+ "/../../" +"correction/POGCorr/POG/MUO/2023_Summer23BPix/muon_scalesmearing.json",
        "Run2022C": this_dir+ "/../../" +"correction/POGCorr/POG/MUO/2022_Summer22/muon_scalesmearing.json",
        "Run2022D": this_dir+ "/../../" +"correction/POGCorr/POG/MUO/2022_Summer22/muon_scalesmearing.json",
        "Run2022E": this_dir+ "/../../" +"correction/POGCorr/POG/MUO/2022_Summer22EE/muon_scalesmearing.json",
        "Run2022F": this_dir+ "/../../" +"correction/POGCorr/POG/MUO/2022_Summer22EE/muon_scalesmearing.jsonn",
        "Run2022G": this_dir+ "/../../" +"correction/POGCorr/POG/MUO/2022_Summer22EE/muon_scalesmearing.json",
        "Run2024C": this_dir+ "/../../" +"correction/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/muon_scalesmearing.json",
        "Run2024D": this_dir+ "/../../" +"correction/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/muon_scalesmearing.json",
        "Run2024E": this_dir+ "/../../" +"correction/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/muon_scalesmearing.json",
        "Run2024F": this_dir+ "/../../" +"correction/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/muon_scalesmearing.json",
        "Run2024G": this_dir+ "/../../" +"correction/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/muon_scalesmearing.json",
        "Run2024H": this_dir+ "/../../" +"correction/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/muon_scalesmearing.json",
        "Run2024I": this_dir+ "/../../" +"correction/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/muon_scalesmearing.json",
    }

    jsonLists = {
        "Run2023C": this_dir+ "/../" +"json/samples/Run2023C.json",
        "Run2023D": this_dir+ "/../" +"json/samples/Run2023D.json",
        "Run2022C": this_dir+ "/../" +"json/samples/Run2022C.json",
        "Run2022D": this_dir+ "/../" +"json/samples/Run2022D.json",
        "Run2022E": this_dir+ "/../" +"json/samples/Run2022E.json",
        "Run2022F": this_dir+ "/../" +"json/samples/Run2022F.json",
        "Run2022G": this_dir+ "/../" +"json/samples/Run2022G.json",

        "Run2024C": this_dir+ "/../" +"json/samples/Run2024C.json",
        "Run2024D": this_dir+ "/../" +"json/samples/Run2024D.json",
        "Run2024E": this_dir+ "/../" +"json/samples/Run2024E.json",
        "Run2024F": this_dir+ "/../" +"json/samples/Run2024F.json",
        "Run2024G": this_dir+ "/../" +"json/samples/Run2024G.json",
        "Run2024H": this_dir+ "/../" +"json/samples/Run2024H.json",
        "Run2024I": this_dir+ "/../" +"json/samples/Run2024I.json",
    }
    datasets = {
        "Run2023C": ["Muon0", "Muon1"], 
        "Run2023D": ["Muon0", "Muon1"],
        "Run2022C": ["Muon"],
        "Run2022D": ["Muon"],
        "Run2022E": ["Muon"],
        "Run2022F": ["Muon"],
        "Run2022G": ["Muon"],

        "Run2024C": ["Muon0", "Muon1"],
        "Run2024D": ["Muon0", "Muon1"],
        "Run2024E": ["Muon0", "Muon1"],
        "Run2024F": ["Muon0", "Muon1"],
        "Run2024F": ["Muon0", "Muon1"],
        "Run2024G": ["Muon0", "Muon1"],
        "Run2024H": ["Muon0", "Muon1"],
        "Run2024I": ["Muon0", "Muon1", "Muon0_add", "Muon1_add"]
    }

    for era in periods:
        print("skim era: " + str(era))
        ROOT.gInterpreter.ProcessLine('load_correction("' + correctionFiles[era] +'");')
        # ROOT.gInterpreter.ProcessLine('bind_correction("' + era + '");')
        # postProcess.correction(jsonLists, era, datasets, RochesterCorrection)
        postProcess.correction(jsonLists, era, datasets[era], RochesterCorrection)

if __name__ == "__main__":
    main()
