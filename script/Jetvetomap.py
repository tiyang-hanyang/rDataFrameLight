import ROOT
import correctionlib
import json
import os

class correctionApplier:
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

    def reweight_and_save(self):
        rdf = ROOT.ROOT.RDataFrame(self.c1)
        print("dataframe created")
        # only record the ZptWgt in the definition
        branchArray = list(rdf.GetColumnNames())
        branchArray.append("PassVetoJet_jetId")
        branchArray.append("PassVetoJet_pt")
        branchArray.append("PassVetoJet_phi")
        branchArray.append("PassVetoJet_eta")
        branchArray.append("PassVetoJet_mass")
        branchArray.append("nPassVetoJet")
        branchArray.append("PassVetoBJet_jetId")
        branchArray.append("PassVetoBJet_pt")
        branchArray.append("PassVetoBJet_phi")
        branchArray.append("PassVetoBJet_eta")
        branchArray.append("PassVetoBJet_mass")
        branchArray.append("nPassVetoBJet")

        # rdf = rdf.Define("PassJetVeto",'corr->evaluate(std::vector<correction::Variable::Type>{"jetvetomap", GoodJet_eta, GoodJet_phi}) == 0')
        # rdf = rdf.Define("PassJetVeto", """
        #     [](const ROOT::VecOps::RVec<float>& eta, const ROOT::VecOps::RVec<float>& phi) -> ROOT::VecOps::RVec<int> {
        #         ROOT::VecOps::RVec<int> pass_flags(eta.size());
        #         for (size_t i = 0; i < eta.size(); ++i) {
        #             pass_flags[i] = (corr->evaluate(std::vector<correction::Variable::Type>{"jetvetomap", eta[i], phi[i]}) == 0);
        #         }
        #         return pass_flags;
        #     }
        # """, ["GoodJet_eta", "GoodJet_phi"])
        rdf = rdf.Define("PassJetVeto", "passJetVetoFunc(GoodJet_eta, GoodJet_phi)")
        rdf = rdf.Define("PassVetoJet_jetId", "GoodJet_jetId[PassJetVeto]")
        rdf = rdf.Define("PassVetoJet_pt", "GoodJet_pt[PassJetVeto]")
        rdf = rdf.Define("PassVetoJet_eta", "GoodJet_eta[PassJetVeto]")
        rdf = rdf.Define("PassVetoJet_phi", "GoodJet_phi[PassJetVeto]")
        rdf = rdf.Define("PassVetoJet_mass", "GoodJet_mass[PassJetVeto]")
        rdf = rdf.Define("nPassVetoJet", "PassVetoJet_jetId.size()")
        # b jet 
        # rdf = rdf.Define("PassJetVetoBJet",'corr->evaluate(std::vector<correction::Variable::Type>{"jetvetomap", BJet_eta, BJet_phi}) == 0')

        rdf = rdf.Define("PassJetVetoBJet", "passJetVetoFunc(BJet_eta, BJet_phi)")
        rdf = rdf.Define("PassVetoBJet_jetId", "BJet_jetId[PassJetVetoBJet]")
        rdf = rdf.Define("PassVetoBJet_pt", "BJet_pt[PassJetVetoBJet]")
        rdf = rdf.Define("PassVetoBJet_eta", "BJet_eta[PassJetVetoBJet]")
        rdf = rdf.Define("PassVetoBJet_phi", "BJet_phi[PassJetVetoBJet]")
        rdf = rdf.Define("PassVetoBJet_mass", "BJet_mass[PassJetVetoBJet]")
        rdf = rdf.Define("nPassVetoBJet", "PassVetoBJet_jetId.size()")
        print("jetveomap DONE")

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
    ROOT.EnableImplicitMT()
    ROOT.gErrorIgnoreLevel = ROOT.kInfo

    periods = [
        # "Run2023C", 
        "Run2023D"
        ]
    jsonLists = {
        "Run2023C": "/home/tiyang/public/rDataFrameLight_git/source/json/samples/MetricMuon/MetricMuonSkimmed_Run2023C.json",
        "Run2023D": "/home/tiyang/public/rDataFrameLight_git/source/json/samples/MetricMuon/MetricMuonSkimmed_Run2023D.json"
    }
    correctionFiles = {
        "Run2023C": "/home/tiyang/public/rDataFrameLight_git/correction/JEC/jsonpog-integration-master/POG/JME/2023_Summer23/jetvetomaps.json",
        "Run2023D": "/home/tiyang/public/rDataFrameLight_git/correction/JEC/jsonpog-integration-master/POG/JME/2023_Summer23BPix/jetvetomaps.json"
    }
    corrName = {
        "Run2023C": "Summer23Prompt23_RunC_V1",
        "Run2023D": "Summer23BPixPrompt23_RunD_V1"
    }
    datasets = [
        # "Muon0",
        # "Muon1",
        # "DY2L",
        # "ttbarDL",
        # "ttbarSL",
        # "WJets",

        # "WW2L2Nu",
        # "WZ2L2Q",
        # "ZZ2L2Nu",
        # "ZZ2L2Q",

        # "WWW",
        # "WWZ",
        # "WZZ",
        # "ZZZ",

        # "TTTT",
        # "TTHBB",
        # "TTHZZ",
        # "TTWH",
        # "TTZH",
        # "TTWW",
        # "TTWZ",
        # "TTZZ",
        "QCD_15_20_mu",
        "QCD_20_30_mu",
        "QCD_30_50_mu",
        "QCD_50_80_mu",
        "QCD_80_120_mu",
        "QCD_120_170_mu",
        "QCD_170_300_mu",
        "QCD_300_470_mu",
        "QCD_470_600_mu",
        "QCD_600_800_mu",
        "QCD_800_1000_mu",
        "QCD_1000_mu",
        "TWm2L",
        "TbarWp2L",

        "TTHH_DL_2B2V",
        "TTHH_DL_2B2Tau",
        "TTHH_SL_2B2V",
        "TTHH_SL_2B2Tau",

        "TTZ_low",
        "TTZ_high",
        "TTW",
        "TTHnonBB",
        "THQ",
        "THW",
        "GGZZ_2mu2e",
        "GGZZ_4mu",
        "GGZZ_4e",
        "GGZZ_2mu2tau",
        "GGZZ_2e2tau",
        "DYG",
        "WG",
        "ZZ4l",
        "WZ3l",
        "TZQB",
    ]
    correctionlib.register_pyroot_binding()

    # using pointer to control the replace of variables (correction configs)
    ROOT.gInterpreter.Declare("""
        #include "correction.h"
        #include <memory>
        std::unique_ptr<correction::CorrectionSet> cset;
        std::shared_ptr<const correction::Correction> corr;

        void load_correction(const std::string& filepath) {
            cset = correction::CorrectionSet::from_file(filepath);
        }

        void bind_correction(const std::string& name) {
            if (!cset) {
                throw std::runtime_error("CorrectionSet not loaded before binding!");
            }
            corr = cset->at(name);
        }
    """)
    # ROOT.gInterpreter.Declare('auto cset = correction::CorrectionSet::from_file("/home/tiyang/public/rDataFrameLight_git/correction/JEC/jsonpog-integration-master/POG/JME/2023_Summer23BPix/jetvetomaps.json");')
    # ROOT.gInterpreter.Declare('auto corr = cset->at("DY_pTll_reweighting");')

    # definition the passJetVetoFunc classifier commonly needed
    # ROOT.gInterpreter.Declare('''
    #     #include <functional>
    #     std::function<ROOT::VecOps::RVec<int>(ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float>)> passJetVetoFunc =
    #     [](ROOT::VecOps::RVec<float> eta, ROOT::VecOps::RVec<float>phi) {
    #         ROOT::VecOps::RVec<int> pass_flags(eta.size());
    #         for (size_t i = 0; i < eta.size(); ++i) {
    #             float safe_phi = std::clamp(phi[i], -3.1415f, 3.1415f);
    #             pass_flags[i] = (corr->evaluate(std::vector<correction::Variable::Type>{"jetvetomap", eta[i], safe_phi}) == 0);
    #         }
    #         return pass_flags;
    #     };
    # ''')
    ROOT.gInterpreter.Declare('''
    #include <functional>
    ROOT::VecOps::RVec<int> passJetVetoFunc(ROOT::VecOps::RVec<float> eta, ROOT::VecOps::RVec<float> phi) {
        ROOT::VecOps::RVec<int> pass_flags(eta.size());
        for (size_t i = 0; i < eta.size(); ++i) {
            float safe_phi = std::clamp(phi[i], -3.1415f, 3.1415f);
            pass_flags[i] = (corr->evaluate(std::vector<correction::Variable::Type>{"jetvetomap", eta[i], safe_phi}) == 0);
        }
        return pass_flags;
    }
    ''')

    for era in periods:
        print("skim era: " + str(era))
        ROOT.gInterpreter.ProcessLine('load_correction("'+ correctionFiles[era] +'");')
        ROOT.gInterpreter.ProcessLine('bind_correction("'+corrName[era]+'");')
        for dataset in datasets:
            print("skim dataset: " + str(dataset))
            eraProcessor = correctionApplier(jsonLists[era], era, dataset)
            if eraProcessor.numberOfFiles == 0:
                continue
            eraProcessor.reweight_and_save()

if __name__ == "__main__":
    main()