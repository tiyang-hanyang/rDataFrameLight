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
        branchArray = list(rdf.GetColumnNames())
        branchArray.append("leadingMuonPt_Rscale")
        branchArray.append("subleadingMuonPt_Rscale")
        branchArray.append("leadingMuonPt_Rcorr")
        branchArray.append("subleadingMuonPt_Rcorr")
        branchArray.append("leadingMuonPt_Rscale_up")
        branchArray.append("leadingMuonPt_Rscale_dn")
        branchArray.append("subleadingMuonPt_Rscale_up")
        branchArray.append("subleadingMuonPt_Rscale_dn")
        branchArray.append("leadingMuonPt_Rcorr_resolup")
        branchArray.append("leadingMuonPt_Rcorr_resoldn")
        branchArray.append("subleadingMuonPt_Rcorr_resolup")
        branchArray.append("subleadingMuonPt_Rcorr_resoldn")

        # scale correction
        rdf = rdf.Define('leadingMuonPt_Rscale','pt_scale(0, leadingMuonPt, leadingMuonEta, leadingMuonPhi, leadingMuonCharge)')
        rdf = rdf.Define('subleadingMuonPt_Rscale','pt_scale(0, subleadingMuonPt, subleadingMuonEta, subleadingMuonPhi, subleadingMuonCharge )')

        # resolution correction
        rdf = rdf.Define('leadingMuonPt_Rcorr','pt_resol(leadingMuonPt_Rscale, leadingMuonEta, float(GoodMuon_nTrackerLayers[0]))')
        rdf = rdf.Define('subleadingMuonPt_Rcorr','pt_resol(subleadingMuonPt_Rscale, subleadingMuonEta, float(GoodMuon_nTrackerLayers[1]))')

        # uncertainties
        rdf = rdf.Define('leadingMuonPt_Rscale_up', 'pt_scale_var(leadingMuonPt_Rcorr, leadingMuonEta, leadingMuonPhi, leadingMuonCharge, "up")')
        rdf = rdf.Define('leadingMuonPt_Rscale_dn', 'pt_scale_var(leadingMuonPt_Rcorr, leadingMuonEta, leadingMuonPhi, leadingMuonCharge, "dn")')
        rdf = rdf.Define('subleadingMuonPt_Rscale_up', 'pt_scale_var(subleadingMuonPt_Rcorr, subleadingMuonEta, subleadingMuonPhi, subleadingMuonCharge, "up")')
        rdf = rdf.Define('subleadingMuonPt_Rscale_dn', 'pt_scale_var(subleadingMuonPt_Rcorr, subleadingMuonEta, subleadingMuonPhi, subleadingMuonCharge, "dn")')

        rdf = rdf.Define("leadingMuonPt_Rcorr_resolup", 'pt_resol_var(leadingMuonPt_Rscale, leadingMuonPt_Rcorr, leadingMuonEta, "up")')
        rdf = rdf.Define("leadingMuonPt_Rcorr_resoldn", 'pt_resol_var(leadingMuonPt_Rscale, leadingMuonPt_Rcorr, leadingMuonEta, "dn")')
        rdf = rdf.Define("subleadingMuonPt_Rcorr_resolup", 'pt_resol_var(subleadingMuonPt_Rscale, subleadingMuonPt_Rcorr, subleadingMuonEta, "up")')
        rdf = rdf.Define("subleadingMuonPt_Rcorr_resoldn", 'pt_resol_var(subleadingMuonPt_Rscale, subleadingMuonPt_Rcorr, subleadingMuonEta, "dn")')

        print("Rochester correction DONE")

        # saving with corresponding names
        print("preparing save options")
        opt = ROOT.RDF.RSnapshotOptions()
        opt.fCompressionLevel = 9

        print("create out folder")
        outDir = "./RochesterCorrection/"+self.era+"/"
        if not os.path.exists(outDir):
            os.makedirs(outDir)
        outputPath = outDir+self.dataset+"_skimmed_Rcorr.root"
        print("save")
        rdf.Snapshot("Events", outputPath, branchArray, opt)
        print("done")

def main():
    postProcess.initializing()
    this_dir = os.path.dirname(os.path.abspath(__file__))
    ROOT.gInterpreter.AddIncludePath(this_dir)
    ROOT.gInterpreter.ProcessLine('#include "MuonScaRe.cc"')

    # specify the periods to run
    periods = [
        "Run3Summer23NanoAODv12", 
        "Run3Summer23BPixNanoAODv12"
    ]

    # correctionlib json files accordingly
    correctionFiles = {
        "Run3Summer23NanoAODv12": "/home/tiyang/public/rDataFrameLight_git/correction/Rochester/muonscarekit-master/corrections/2023_Summer23.json",
        "Run3Summer23BPixNanoAODv12": "/home/tiyang/public/rDataFrameLight_git/correction/Rochester/muonscarekit-master/corrections/2023_Summer23BPix.json"
    }

    from MetricSkimmedFiles import RochesterFileLists as jsonLists
    from MetricSkimmedFiles import datasets

    for era in periods:
        print("skim era: " + str(era))
        ROOT.gInterpreter.ProcessLine('load_correction("' + correctionFiles[era] +'");')
        # ROOT.gInterpreter.ProcessLine('bind_correction("' + era + '");')
        postProcess.correction(jsonLists, era, datasets[era], RochesterCorrection)

if __name__ == "__main__":
    main()