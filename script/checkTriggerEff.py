## using only dimuon trigger
import ROOT 
import os
from array import array

ROOT.EnableImplicitMT()

def getRDF(type):
    rep = {
        "DL": "/data2/common/NanoAOD/private_v15/mc/RunIII2024Summer24NanoAODv15/TTHH_DL_LO_2B2W_251024_002450/",
        "SL": "/data2/common/NanoAOD/private_v15/mc/RunIII2024Summer24NanoAODv15/TTHH_SL_LO_2B2W_251024_082942/",
    }
    allFiles = [rep[type]+"0000/"+fName for fName in os.listdir(rep[type]+"0000/")] + [rep[type]+"0001/"+fName for fName in os.listdir(rep[type]+"0001/")]
    allFiles = list([fName for fName in allFiles if "log" not in fName])
    rdf = ROOT.RDataFrame("Events", allFiles)
    rdf = rdf.Define("GoodMuon", "(Muon_pt>15)*(abs(Muon_eta)<2.4) && (abs(Muon_dxy)<0.05) && (abs(Muon_dz)<0.1) && (abs(Muon_sip3d)<8) && Muon_mediumId && (Muon_miniPFRelIso_all<0.4) && (Muon_jetDF<0.2480) && (Muon_promptMVA > 0.64)")
    rdf = rdf.Define("nGoodMuon", "Nonzero(GoodMuon).size()").Filter("nGoodMuon>1")
    rdf = rdf.Define("leadingMuon_pt", "Muon_pt[Nonzero(GoodMuon)[0]]")
    return rdf

def plotHistLinear(rdf, type, trigger, triggerText, outName):
    pTModel = ROOT.RDF.TH1DModel("hPassTr", "", 17, 15., 100.),
    passTr = rdf.Filter(trigger).Histo1D(pTModel, "leadingMuon_pt")
    pTModel2 = ROOT.RDF.TH1DModel("hAll", "",  17, 15., 100.)
    histAll = rdf.Histo1D(pTModel2, "leadingMuon_pt")

    # ratio
    ROOT.RDF.RunGraphs([passTr, histAll])
    hnum = passTr.GetValue()
    hden = histAll.GetValue()
    triggerEff = hnum.Clone("eff")
    triggerEff.Divide(hnum, hden, 1.0, 1.0, "B")
    triggerEff.SetStats(0)
    triggerEff.SetTitle("")
    triggerEff.GetYaxis().SetTitle("Trigger eff")
    triggerEff.GetXaxis().SetTitle("leading Muon p_{T} [GeV]")
    triggerEff.SetLineWidth(2)
    triggerEff.SetMaximum(1.0)
    triggerEff.SetMinimum(0.5)
    # bin values
    eff_values = []
    eff_errors = []
    for i in range(1, triggerEff.GetNbinsX() + 1):
        eff = triggerEff.GetBinContent(i)
        err = triggerEff.GetBinError(i)
        eff_values.append(eff)
        eff_errors.append(err)

    # save (with writing)
    c1 = ROOT.TCanvas("c1", "c1", 600, 700)
    c1.SetLeftMargin(0.14)
    triggerEff.Draw("HIST E")
    latex = ROOT.TLatex()
    latex.SetTextSize(0.03)
    latex.SetTextFont(42)
    latex.SetTextAlign(12)
    latex.SetNDC(True)
    latex.DrawLatex(0.2, 0.85, type)
    startingY=0.8
    for i in triggerText:
        latex.DrawLatex(0.2, startingY, i)
        startingY -= 0.05
    c1.Update()
    c1.SaveAs(outName)


def plotHists(rdf, type, trigger, triggerText, outName):
    pTModel = ROOT.RDF.TH1DModel("hPassTr", "", 10, array('d',[15,16,17,18,19,20,21,25,30,50,100]))
    passTr = rdf.Filter(trigger).Histo1D(pTModel, "leadingMuon_pt")
    pTModel2 = ROOT.RDF.TH1DModel("hAll", "", 10, array('d',[15,16,17,18,19,20,21,25,30,50,100]))
    histAll = rdf.Histo1D(pTModel2, "leadingMuon_pt")

    # ratio
    ROOT.RDF.RunGraphs([passTr, histAll])
    hnum = passTr.GetValue()
    hden = histAll.GetValue()
    triggerEff = hnum.Clone("eff")
    triggerEff.Divide(hnum, hden, 1.0, 1.0, "B")

    # bin hist
    labelBinHist = ROOT.TH1D("eff_bin", "eff_bin", 10, 0.0, 10.0)
    labelBinHist.SetStats(0)
    labelBinHist.SetTitle("")
    labelBinHist.GetYaxis().SetTitle("Trigger eff")
    labelBinHist.GetXaxis().SetTitle("leading Muon p_{T} [GeV]")
    labelBinHist.GetXaxis().SetBinLabel(1, "15-16")
    labelBinHist.GetXaxis().SetBinLabel(2, "16-17")
    labelBinHist.GetXaxis().SetBinLabel(3, "17-18")
    labelBinHist.GetXaxis().SetBinLabel(4, "18-19")
    labelBinHist.GetXaxis().SetBinLabel(5, "19-20")
    labelBinHist.GetXaxis().SetBinLabel(6, "20-21")
    labelBinHist.GetXaxis().SetBinLabel(7, "21-25")
    labelBinHist.GetXaxis().SetBinLabel(8, "25-30")
    labelBinHist.GetXaxis().SetBinLabel(9, "30-50")
    labelBinHist.GetXaxis().SetBinLabel(10, "50-100")
    labelBinHist.SetLineWidth(2)
    labelBinHist.SetMaximum(1.0)
    labelBinHist.SetMinimum(0.2)
    # bin values
    eff_values = []
    eff_errors = []
    for i in range(1, triggerEff.GetNbinsX() + 1):
        eff = triggerEff.GetBinContent(i)
        labelBinHist.SetBinContent(i, eff)
        err = triggerEff.GetBinError(i)
        labelBinHist.SetBinError(i, err)
        eff_values.append(eff)
        eff_errors.append(err)

    # save (with writing)
    c1 = ROOT.TCanvas("c1", "c1", 600, 700)
    labelBinHist.Draw("HIST E")
    latex = ROOT.TLatex()
    latex.SetTextAlign(22) 
    latex.SetTextSize(0.025)
    latex.SetTextFont(42)
    latex.SetNDC(False)
    for i, (eff, err) in enumerate(zip(eff_values, eff_errors), start=1):
        x = labelBinHist.GetBinCenter(i)
        y = labelBinHist.GetBinContent(i)
        y_text = y - 0.1
        if y_text < 0.2:
            y_text = 0.3
        latex.DrawLatex(x, y_text, f"{eff:.2f}")
        latex.DrawLatex(x, y_text-0.05, f"#pm {err:.2f}")
    latex.SetTextAlign(32)
    latex.SetNDC(True)
    startingY=0.5
    for i in triggerText:
        latex.DrawLatex(0.85, startingY, i)
        startingY -= 0.05
    latex.SetTextAlign(12)
    latex.DrawLatex(0.1, 0.85, type)
    c1.Update()
    c1.SaveAs(outName)

if __name__ == "__main__":
    rdf = getRDF("SL")
    # plotHists(rdf, "TTHH SL", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"],"TTHH_SL_triggerEff_dimuonRec.png")
    # plotHists(rdf, "TTHH SL", "(HLT_IsoMu24 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8)", ["HLT_IsoMu24", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"], "TTHH_SL_triggerEff_singleDimuonRec.png")
    # plotHists(rdf, "TTHH SL", "(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_IsoMu24 || HLT_IsoMu24_eta2p1 || HLT_IsoMu27 || HLT_Mu50 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8)", ["full single, di-muon triggers"],"TTHH_SL_triggerEff_singleDimuonAll.png")
    # plotHists(rdf, "TTHH SL", "(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_IsoMu24 || HLT_IsoMu24_eta2p1 || HLT_IsoMu27 || HLT_Mu50 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8 || HLT_TripleMu_5_3_3_Mass3p8_DZ || HLT_TripleMu_10_5_5_DZ || HLT_TripleMu_12_10_5)", ["full single, di, tri-muon triggers"],"TTHH_SL_triggerEff_full.png")

    plotHistLinear(rdf, "TTHH SL", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"],"TTHH_SL_triggerEff_dimuonRec.png")
    plotHistLinear(rdf, "TTHH SL", "(HLT_IsoMu24 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8)", ["HLT_IsoMu24", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"], "TTHH_SL_triggerEff_singleDimuonRec.png")
    plotHistLinear(rdf, "TTHH SL", "(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_IsoMu24 || HLT_IsoMu24_eta2p1 || HLT_IsoMu27 || HLT_Mu50 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8)", ["Single and di-muon triggers"],"TTHH_SL_triggerEff_singleDimuonAll.png")
    plotHistLinear(rdf, "TTHH SL", "(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_IsoMu24 || HLT_IsoMu24_eta2p1 || HLT_IsoMu27 || HLT_Mu50 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8 || HLT_TripleMu_5_3_3_Mass3p8_DZ || HLT_TripleMu_10_5_5_DZ || HLT_TripleMu_12_10_5)", ["Full trigger menu"],"TTHH_SL_triggerEff_full.png")

    rdf = getRDF("DL")
    # plotHists(rdf, "TTHH DL", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"],"TTHH_DL_triggerEff_dimuonRec.png")
    # plotHists(rdf, "TTHH DL", "(HLT_IsoMu24 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8)", ["HLT_IsoMu24", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"], "TTHH_DL_triggerEff_singleDimuonRec.png")
    # plotHists(rdf, "TTHH DL", "(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_IsoMu24 || HLT_IsoMu24_eta2p1 || HLT_IsoMu27 || HLT_Mu50 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8)", ["full single, di-muon triggers"],"TTHH_DL_triggerEff_singleDimuonAll.png")
    # plotHists(rdf, "TTHH DL", "(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_IsoMu24 || HLT_IsoMu24_eta2p1 || HLT_IsoMu27 || HLT_Mu50 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8 || HLT_TripleMu_5_3_3_Mass3p8_DZ || HLT_TripleMu_10_5_5_DZ || HLT_TripleMu_12_10_5)", ["full single, di, tri-muon triggers"],"TTHH_DL_triggerEff_full.png")

    plotHistLinear(rdf, "TTHH DL", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"],"TTHH_DL_triggerEff_dimuonRec.png")
    plotHistLinear(rdf, "TTHH DL", "(HLT_IsoMu24 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8)", ["HLT_IsoMu24", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"], "TTHH_DL_triggerEff_singleDimuonRec.png")
    plotHistLinear(rdf, "TTHH DL", "(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_IsoMu24 || HLT_IsoMu24_eta2p1 || HLT_IsoMu27 || HLT_Mu50 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8)", ["Single and di-muon triggers"],"TTHH_DL_triggerEff_singleDimuonAll.png")
    plotHistLinear(rdf, "TTHH DL", "(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_IsoMu24 || HLT_IsoMu24_eta2p1 || HLT_IsoMu27 || HLT_Mu50 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8 || HLT_TripleMu_5_3_3_Mass3p8_DZ || HLT_TripleMu_10_5_5_DZ || HLT_TripleMu_12_10_5)", ["Full trigger menu"],"TTHH_DL_triggerEff_full.png")