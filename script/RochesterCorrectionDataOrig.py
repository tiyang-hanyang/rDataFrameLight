import postProcess
import ROOT
import correctionlib
import os

# redefine the branch operation method and directory to save
class RochesterCorrection(postProcess.correctionApplier):
    def reweight_and_save(self):
        # preparing the computation of JetVetoMap
        rdf = ROOT.ROOT.RDataFrame(self.c1)
        print("dataframe created")
        #branchArray = list(rdf.GetColumnNames())
        branchArray = [
            "run",
            "luminosityBlock",
            "event",

            "PV_npvsGood",

            "nJet",
            "Jet_jetId",
            "Jet_area",
            "Jet_pt",
            "Jet_eta",
            "Jet_phi",
            "Jet_mass",
            "Jet_btagPNetB",

            "nMuon",
            "Muon_pt",
            "Muon_ptErr",
            "Muon_eta",
            "Muon_phi",
            "Muon_mass",
            "Muon_charge",
            "Muon_miniPFRelIso_all",
            "Muon_pfRelIso04_all",
            "Muon_tightId",
            "Muon_dxy",
            "Muon_dxyErr",
            "Muon_dz",
            "Muon_dzErr",
            "Muon_nTrackerLayers",

            "MET_phi",
            "MET_pt",
            "MET_sumEt",

            "Flag_goodVertices",
            "Flag_globalSuperTightHalo2016Filter",
            "Flag_EcalDeadCellTriggerPrimitiveFilter",
            "Flag_BadPFMuonFilter",
            "Flag_eeBadScFilter",
            "Flag_BadPFMuonDzFilter",
            "Flag_hfNoisyHitsFilter",

            "HLT_IsoMu20",
            "HLT_IsoMu24",
            "HLT_IsoMu24_eta2p1",
            "HLT_IsoMu27",
            "HLT_Mu50",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
            "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8",
            "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8",
            "HLT_TripleMu_5_3_3_Mass3p8_DZ",
            "HLT_TripleMu_10_5_5_DZ",
            "HLT_TripleMu_12_10_5",

            "HLT_CascadeMu100",
            "HLT_HighPtTkMu100"
        ]
        branchArray.append("Muon_pt_Rscale")

        # data case
        rdf = rdf.Define('Muon_pt_Rscale', 'pt_scale(1, Muon_pt, Muon_eta, Muon_phi, Muon_charge)')

        print("Rochester correction DONE")

        # saving with corresponding names
        print("preparing save options")
        opt = ROOT.RDF.RSnapshotOptions()
        opt.fCompressionLevel = 9

        print("create out folder")
        outDir = "./RochesterCorrection/data/"+self.era+"/"+self.dataset+"/"
        if not os.path.exists(outDir):
            os.makedirs(outDir)
        outputPath = outDir+self.dataset+"_Rcorr.root"
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
        "Run2023C", 
        "Run2023D"
    ]

    # correctionlib json files accordingly
    correctionFiles = {
        "Run2023C": "/home/tiyang/public/rDataFrameLight_git/correction/Rochester/muonscarekit-master/corrections/2023_Summer23.json",
        "Run2023D": "/home/tiyang/public/rDataFrameLight_git/correction/Rochester/muonscarekit-master/corrections/2023_Summer23BPix.json"
    }

    jsonLists = {
        "Run2023C": "/home/tiyang/public/rDataFrameLight_git/source/json/samples/Run2023C.json",
        "Run2023D": "/home/tiyang/public/rDataFrameLight_git/source/json/samples/Run2023D.json"
    }
    datasets = ["Muon0", "Muon1"]

    for era in periods:
        print("skim era: " + str(era))
        ROOT.gInterpreter.ProcessLine('load_correction("' + correctionFiles[era] +'");')
        # ROOT.gInterpreter.ProcessLine('bind_correction("' + era + '");')
        postProcess.correction(jsonLists, era, datasets, RochesterCorrection)

if __name__ == "__main__":
    main()