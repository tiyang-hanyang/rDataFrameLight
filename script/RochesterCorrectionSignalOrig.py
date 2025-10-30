import postProcess
import ROOT
import os
import json

# redefine the branch operation method and directory to save
class RochesterCorrection(postProcess.correctionApplier):
    def reweight_and_save(self):
        for filePath in self.filePaths:
            # first check the existence of output!
            print("create out folder")
            outDir = "/data2/common/skimmed_NanoAOD/skim_1024_RochesterCorr_TTHH/mc/"+self.era+"/"+self.dataset+"/"
            if not os.path.exists(outDir):
                os.makedirs(outDir)

            newFilePath = "-".join(filePath.split("/")[-2:])[:-5]+"_Rcorr.root"
            outputPath = outDir + newFilePath

            if os.path.isfile(outputPath):
                print(outputPath, "exist, skip")
                continue

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
            branchArray.append("Muon_pt_Rscale")
            branchArray.append("Muon_pt_Rcorr")
            branchArray.append("Muon_pt_Rscale_up")
            branchArray.append("Muon_pt_Rscale_dn")
            branchArray.append("Muon_pt_Rcorr_resolup")
            branchArray.append("Muon_pt_Rcorr_resoldn")

            # scale correction
            rdf = rdf.Define('Muon_pt_Rscale', 'pt_scale(0, Muon_pt, Muon_eta, Muon_phi, Muon_charge)')
            # resolution correction
            rdf = rdf.Define('Muon_pt_Rcorr','pt_resol(Muon_pt_Rscale, Muon_eta, Muon_phi, Muon_nTrackerLayers, event, luminosityBlock)')
            # uncertainties
            rdf = rdf.Define('Muon_pt_Rscale_up', 'pt_scale_var(Muon_pt_Rcorr, Muon_eta, Muon_phi, Muon_charge, "up")')
            rdf = rdf.Define('Muon_pt_Rscale_dn', 'pt_scale_var(Muon_pt_Rcorr, Muon_eta, Muon_phi, Muon_charge, "dn")')
            rdf = rdf.Define("Muon_pt_Rcorr_resolup", 'pt_resol_var(Muon_pt_Rscale, Muon_pt_Rcorr, Muon_eta, "up")')
            rdf = rdf.Define("Muon_pt_Rcorr_resoldn", 'pt_resol_var(Muon_pt_Rscale, Muon_pt_Rcorr, Muon_eta, "dn")')

            # for MC, I also need weights
            sum_weights   = rdf.Sum("genWeight").GetValue()

            # saving with corresponding names
            print("preparing save options")
            opt = ROOT.RDF.RSnapshotOptions()
            opt.fCompressionLevel = 9

            print("save")
            rdf.Snapshot("Events", outputPath, branchArray, opt)
            print("done")

            # save all weights
            fRC = ROOT.TFile.Open(outputPath, "UPDATE")
            h_sumw = ROOT.TH1D("genWeightSum", "sum of genWeight (this file)", 1, 0.0, 1.0)
            h_sumw.SetBinContent(1, sum_weights)
            h_sumw.GetYaxis().SetTitle("sum(genWeight)")
            h_sumw.Write("", ROOT.TObject.kOverwrite)
            fRC.Close()

def main():
    postProcess.initializing()
    this_dir = os.path.dirname(os.path.abspath(__file__))
    ROOT.gInterpreter.AddIncludePath(this_dir)
    ROOT.gInterpreter.ProcessLine('#include "MuonScaReOrig.cc"')

    # specify the periods to run
    periods = [
        #"Run3Summer23NanoAODv12", 
        #"Run3Summer23BPixNanoAODv12",
        #"Run3Summer22NanoAODv12",
        #"Run3Summer22EENanoAODv12",
        "RunIII2024Summer24NanoAODv15"
    ]

    # correctionlib json files accordingly
    correctionFiles = {
        "Run3Summer23NanoAODv12": this_dir+ "/../../" +"correction/POGCorr/POG/MUO/2023_Summer23/muon_scalesmearing.json",
        "Run3Summer23BPixNanoAODv12": this_dir+ "/../../" +"correction/POGCorr/POG/MUO/2023_Summer23BPix/muon_scalesmearing.json",
        "Run3Summer22NanoAODv12": this_dir+ "/../../" +"correction/POGCorr/POG/MUO/2022_Summer22/muon_scalesmearing.json",
        "Run3Summer22EENanoAODv12": this_dir+ "/../../" +"correction/POGCorr/POG/MUO/2022_Summer22EE/muon_scalesmearing.json",
        "RunIII2024Summer24NanoAODv15": this_dir+ "/../../" +"correction/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/muon_scalesmearing.json",
    }

    jsonLists = {
        "Run3Summer23NanoAODv12": this_dir+ "/../" +"json/samples/Run3Summer23NanoAODv12.json",
        "Run3Summer23BPixNanoAODv12": this_dir+ "/../" +"json/samples/Run3Summer23BPixNanoAODv12.json",
        "Run3Summer22NanoAODv12": this_dir+ "/../" +"json/samples/Run3Summer22NanoAODv12.json",
        "Run3Summer22EENanoAODv12": this_dir+ "/../" +"json/samples/Run3Summer22EENanoAODv12.json",
        "RunIII2024Summer24NanoAODv15": this_dir+ "/../" +"json/samples/Signal_RunIII2024Summer24NanoAODv15.json",
    }

    # from MetricSkimmedFiles import datasets
    #datasets = ["DY2L", "DY2L_low"]
    datasets = [
        "TTHH_SL_2B2W_batch1",
        "TTHH_SL_2B2W_batch2",
    ]

    for era in periods:
        print("skim era: " + str(era))
        ROOT.gInterpreter.ProcessLine('load_correction("' + correctionFiles[era] +'");')
        # ROOT.gInterpreter.ProcessLine('bind_correction("' + era + '");')
        #postProcess.correction(jsonLists, era, datasets[era], RochesterCorrection)
        postProcess.correction(jsonLists, era, datasets, RochesterCorrection)

if __name__ == "__main__":
    main()
