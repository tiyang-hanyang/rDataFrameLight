import ROOT
import correctionlib

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
        branchArray.append("leadingMuonPt_Rcorr")
        branchArray.append("subleadingMuonPt_Rcorr")

        # data case
        rdf = rdf.Define('leadingMuonPt_Rcorr','pt_scale(1, leadingMuonPt, leadingMuonEta, leadingMuonPhi, GoodMuon_charge[0])')
        rdf = rdf.Define('subleadingMuonPt_Rcorr','pt_scale(1, subleadingMuonPt, subleadingMuonEta, subleadingMuonPhi, GoodMuon_charge[1])')

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
        "Run2023C",
        "Run2023D"
    ]

    # correctionlib json files accordingly
    correctionFiles = {
        "Run2023C": "/home/tiyang/public/rDataFrameLight_git/correction/Rochester/muonscarekit-master/corrections/2023_Summer23.json",
        "Run2023D": "/home/tiyang/public/rDataFrameLight_git/correction/Rochester/muonscarekit-master/corrections/2023_Summer23BPix.json"
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
















# add for tthh
mcList = {
    "Run2023C": "",
    "Run2023D": ""
}
dataList = {
    "Run2023C": "",
    "Run2023D": ""
}
# jsonLists = {
#     "Run2023C": "/home/tiyang/public/rDataFrameLight_git/source/json/samples/MetricMuon/MetricMuonSkimmed_Run2023C.json",
#     "Run2023D": "/home/tiyang/public/rDataFrameLight_git/source/json/samples/MetricMuon/MetricMuonSkimmed_Run2023D.json"
# }
correctionFiles = {
    "Run2023C": "",
    "Run2023D": ""
}

era = "Run2023C"

# done

correctionlib.register_pyroot_binding()
ROOT.EnableImplicitMT()

# tuning importing files
ROOT.gROOT.ProcessLine(
    'auto cset = correction::CorrectionSet::from_file("'+correctionFiles[era]+'");'
)
ROOT.gROOT.ProcessLine('#include "MuonScaRe.cc"')

df_data = ROOT.RDataFrame("Events", dataList[era])
df_mc = ROOT.RDataFrame("Events", mcList[era])

df_data = df_data.Filter("nGoodMuon >=2")
df_mc = df_mc.Filter("nGoodMuon >=2")
# done

# Data apply scale correction
df_data = df_data.Define(
    'pt_1_corr',
    'pt_scale(1, pt_1, eta_1, phi_1, charge_1)'
)

# MC apply scale shift and resolution smearing
df_mc = df_mc.Define(
    'pt_1_scale_corr',
    'pt_scale(0, pt_1, eta_1, phi_1, charge_1)'
)
df_mc = df_mc.Define(
    "pt_1_corr",
    'pt_resol(pt_1_scale_corr, eta_1, float(nTrkLayers_1))'
)

# MC evaluate scale uncertainty
df_mc = df_mc.Define(
    'pt_1_scale_corr_up',
    'pt_scale_var(pt_1_corr, eta_1, phi_1, charge_1, "up")'
)
df_mc = df_mc.Define(
    'pt_1_scale_corr_dn',
    'pt_scale_var(pt_1_corr, eta_1, phi_1, charge_1, "dn")'
)

# MC evaluate resolution uncertainty
df_mc = df_mc.Define(
    "pt_1_corr_resolup",
    'pt_resol_var(pt_1_scale_corr, pt_1_corr, eta_1, "up")'
)
df_mc = df_mc.Define(
    "pt_1_corr_resoldn",
    'pt_resol_var(pt_1_scale_corr, pt_1_corr, eta_1, "dn")'
)

# save results