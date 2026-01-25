## using only dimuon trigger
import ROOT 
import os
from array import array
ROOT.EnableImplicitMT()
rep = "/data2/common/NanoAOD/private_v15/mc/RunIII2024Summer24NanoAODv15/TTHH_SL_LO_2B2W_251024_082942/0000/"
allFiles = ["/data2/common/NanoAOD/private_v15/mc/RunIII2024Summer24NanoAODv15/TTHH_SL_LO_2B2W_251024_082942/0000/"+fName for fName in os.listdir("/data2/common/NanoAOD/private_v15/mc/RunIII2024Summer24NanoAODv15/TTHH_SL_LO_2B2W_251024_082942/0000/")] + ["/data2/common/NanoAOD/private_v15/mc/RunIII2024Summer24NanoAODv15/TTHH_SL_LO_2B2W_251024_082942/0001/"+fName for fName in os.listdir("/data2/common/NanoAOD/private_v15/mc/RunIII2024Summer24NanoAODv15/TTHH_SL_LO_2B2W_251024_082942/0001/")]
allFiles = list([fName for fName in allFiles if "log" not in fName])
rdf = ROOT.RDataFrame("Events", allFiles)
rdf = rdf.Define("GoodMuon", "(Muon_pt>15) * (abs(Muon_eta)<2.4) && (abs(Muon_dxy)<0.05) && (abs(Muon_dz)<0.1) && (abs(Muon_sip3d)<8) && Muon_mediumId && (Muon_miniPFRelIso_all<0.4) && (Muon_jetDF<0.2480) && (Muon_promptMVA > 0.64)")
rdf = rdf.Define("nGoodMuon", "Nonzero(GoodMuon).size()").Filter("nGoodMuon>1")
rdf = rdf.Define("leadingMuon_pt", "Muon_pt[Nonzero(GoodMuon)[0]]")


# # if with only two cases
# rdf = ROOT.RDataFrame("Events", allFiles)
# rdf = rdf.Define("GoodMuon", "(Muon_pt>15) && (abs(Muon_eta)<2.4) && (abs(Muon_dxy)<0.05) && (abs(Muon_dz)<0.1) && (abs(Muon_sip3d)<8) && Muon_mediumId && (Muon_miniPFRelIso_all<0.4) && (Muon_jetDF<0.2480) && (Muon_promptMVA > 0.64)")
# rdf = rdf.Define("nGoodMuon", "Nonzero(GoodMuon).size()").Filter("nGoodMuon==2")
# rdf = rdf.Define("leadingMuon_pt", "Muon_pt[Nonzero(GoodMuon)[0]]")

# get hists
pTModel = ROOT.RDF.TH1DModel("hPassTr", "", 8, array('d',[15,17,19,21,25,30,50,100,1000]))
passTr = rdf.Filter("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8").Histo1D(pTModel, "leadingMuon_pt")
pTModel2 = ROOT.RDF.TH1DModel("hAll", "", 8, array('d',[15,17,19,21,25,30,50,100,1000]))
histAll = rdf.Histo1D(pTModel2, "leadingMuon_pt")
# ratio
ROOT.RDF.RunGraphs([passTr, histAll])
hnum = passTr.GetValue()
hden = histAll.GetValue()
triggerEff = hnum.Clone("eff")
triggerEff.Divide(hden)
triggerEff.GetYaxis().SetTitle("Trigger eff")
triggerEff.GetXaxis().SetTitle("leading Muon p_{T} [GeV]")
triggerEff.SetStats(0)
# triggerEff.GetXaxis().SetBinLabel(1, "15-17")
# triggerEff.GetXaxis().SetBinLabel(2, "17-19")
# triggerEff.GetXaxis().SetBinLabel(3, "19-21")
# triggerEff.GetXaxis().SetBinLabel(4, "21-25")
# triggerEff.GetXaxis().SetBinLabel(5, "25-30")
# triggerEff.GetXaxis().SetBinLabel(6, "30-50")
# triggerEff.GetXaxis().SetBinLabel(7, "50-100")
# triggerEff.GetXaxis().SetBinLabel(8, "100-1000")
triggerEff.SetLineWidth(2)
triggerEff.SetMaximum(1.0)
triggerEff.SetMinimum(0.2)
# values
# eff_values = []
# for i in range(1, hnum.GetNbinsX() + 1):
#     n_all  = hden.GetBinContent(i)
#     n_pass = hnum.GetBinContent(i)
#     if n_all > 0:
#         eff = n_pass / n_all
#     else:
#         eff = 0.0
#     eff_values.append(eff)
# save (with writing)
c1 = ROOT.TCanvas("c1", "c1", 600, 700)
c1.SetLogx()
triggerEff.Draw("HIST")
latex = ROOT.TLatex()
latex.SetTextAlign(22) 
latex.SetTextSize(0.03)
latex.SetTextFont(42)
latex.SetNDC(True)
latex.DrawLatex(0.5, 0.8, "only HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8")
latex.SetNDC(False)
# # writing eff
# # latex = ROOT.TLatex()
# # latex.SetTextAlign(22) 
# # latex.SetTextSize(0.03)
# # latex.SetTextFont(42)
# for i, eff in enumerate(eff_values, start=1):
#     x = triggerEff.GetBinCenter(i)
#     y = triggerEff.GetBinContent(i)
#     y_text = 0.9 * y
#     y_text = min(0.98, 0.9*y)
#     latex.DrawLatex(x, y_text, f"{eff:.2f}")
c1.Update()
c1.SaveAs("TTHH_SL_atleast_2mu_triggerEff.png")




# get hists
# Iso Mu 20 is a prescaled trigger! should remove! Looking at prescale = 210
# Mu50
pTModel = ROOT.RDF.TH1DModel("hPassTr", "", 8, array('d',[15,17,19,21,25,30,50,100,1000]))
passTr = rdf.Filter("(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_IsoMu24 || HLT_IsoMu24_eta2p1 || HLT_IsoMu27 || HLT_Mu50 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8 || HLT_TripleMu_5_3_3_Mass3p8_DZ || HLT_TripleMu_10_5_5_DZ || HLT_TripleMu_12_10_5)").Histo1D(pTModel, "leadingMuon_pt")
pTModel2 = ROOT.RDF.TH1DModel("hAll", "", 8, array('d',[15,17,19,21,25,30,50,100,1000]))
histAll = rdf.Histo1D(pTModel2, "leadingMuon_pt")
# ratio
ROOT.RDF.RunGraphs([passTr, histAll])
hnum = passTr.GetValue()
hden = histAll.GetValue()
triggerEff = hnum.Clone("eff")
triggerEff.Divide(hden)
triggerEff.GetYaxis().SetTitle("Trigger eff")
triggerEff.GetXaxis().SetTitle("leading Muon p_{T} [GeV]")
triggerEff.SetStats(0)
# triggerEff.GetXaxis().SetBinLabel(1, "15-17")
# triggerEff.GetXaxis().SetBinLabel(2, "17-19")
# triggerEff.GetXaxis().SetBinLabel(3, "19-21")
# triggerEff.GetXaxis().SetBinLabel(4, "21-25")
# triggerEff.GetXaxis().SetBinLabel(5, "25-30")
# triggerEff.GetXaxis().SetBinLabel(6, "30-50")
# triggerEff.GetXaxis().SetBinLabel(7, "50-100")
# triggerEff.GetXaxis().SetBinLabel(8, "100-1000")
triggerEff.SetLineWidth(2)
triggerEff.SetMaximum(1.0)
triggerEff.SetMinimum(0.2)
c1 = ROOT.TCanvas("c1", "c1", 600, 700)
c1.SetLogx()
triggerEff.Draw("HIST")
latex = ROOT.TLatex()
latex.SetTextAlign(22) 
latex.SetTextSize(0.03)
latex.SetTextFont(42)
latex.SetNDC(True)
latex.DrawLatex(0.5, 0.8, "Full trigger menu")
latex.SetNDC(False)
c1.Update()
c1.SaveAs("TTHH_SL_atleast_2mu_allTriggerEff.png")

# get hists
# pTModel = ROOT.RDF.TH1DModel("hPassTr", "", 8, array('d',[15,17,19,21,25,30,50,100,1000]))
pTModel = ROOT.RDF.TH1DModel("hPassTr", "", 15, array('d',[15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]))
passTr = rdf.Filter("(HLT_IsoMu24 || HLT_IsoMu24_eta2p1 || HLT_IsoMu27 || HLT_Mu50 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8)").Histo1D(pTModel, "leadingMuon_pt")
# pTModel2 = ROOT.RDF.TH1DModel("hPassTr", "", 8, array('d',[15,17,19,21,25,30,50,100,1000]))
pTModel2 = ROOT.RDF.TH1DModel("hPassTr", "", 15, array('d',[15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]))
histAll = rdf.Histo1D(pTModel2, "leadingMuon_pt")
# ratio
ROOT.RDF.RunGraphs([passTr, histAll])
hnum = passTr.GetValue()
hden = histAll.GetValue()
triggerEff = hnum.Clone("eff")
triggerEff.Divide(hden)
triggerEff.GetYaxis().SetTitle("Trigger eff")
triggerEff.GetXaxis().SetTitle("leading Muon p_{T} [GeV]")
triggerEff.SetStats(0)
# triggerEff.GetXaxis().SetBinLabel(1, "15-17")
# triggerEff.GetXaxis().SetBinLabel(2, "17-19")
# triggerEff.GetXaxis().SetBinLabel(3, "19-21")
# triggerEff.GetXaxis().SetBinLabel(4, "21-25")
# triggerEff.GetXaxis().SetBinLabel(5, "25-30")
# triggerEff.GetXaxis().SetBinLabel(6, "30-50")
# triggerEff.GetXaxis().SetBinLabel(7, "50-100")
# triggerEff.GetXaxis().SetBinLabel(8, "100-1000")
triggerEff.SetLineWidth(2)
triggerEff.SetMaximum(1.0)
triggerEff.SetMinimum(0.2)
c1 = ROOT.TCanvas("c1", "c1", 600, 700)
#c1.SetLogx()
triggerEff.Draw("HIST")
latex = ROOT.TLatex()
latex.SetTextAlign(22) 
latex.SetTextSize(0.03)
latex.SetTextFont(42)
latex.SetNDC(True)
latex.DrawLatex(0.5, 0.8, "Full single and di-muon triggers")
latex.SetNDC(False)
c1.Update()
#c1.SaveAs("TTHH_SL_atleast_2mu_singleAndDimuonTriggerEff.png")
c1.SaveAs("TTHH_SL_atleast_2mu_singleAndDimuonTriggerEff_zoomin.png")






# get hists
pTModel = ROOT.RDF.TH1DModel("hPassTr", "", 8, array('d',[15,17,19,21,25,30,50,100,1000]))
passTr = rdf.Filter("(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_IsoMu24)").Histo1D(pTModel, "leadingMuon_pt")
pTModel2 = ROOT.RDF.TH1DModel("hAll", "", 8, array('d',[15,17,19,21,25,30,50,100,1000]))
histAll = rdf.Histo1D(pTModel2, "leadingMuon_pt")
# ratio
ROOT.RDF.RunGraphs([passTr, histAll])
hnum = passTr.GetValue()
hden = histAll.GetValue()
triggerEff = hnum.Clone("eff")
triggerEff.Divide(hden)
triggerEff.GetYaxis().SetTitle("Trigger eff")
triggerEff.GetXaxis().SetTitle("leading Muon p_{T} [GeV]")
triggerEff.SetStats(0)
# triggerEff.GetXaxis().SetBinLabel(1, "15-17")
# triggerEff.GetXaxis().SetBinLabel(2, "17-19")
# triggerEff.GetXaxis().SetBinLabel(3, "19-21")
# triggerEff.GetXaxis().SetBinLabel(4, "21-25")
# triggerEff.GetXaxis().SetBinLabel(5, "25-30")
# triggerEff.GetXaxis().SetBinLabel(6, "30-50")
# triggerEff.GetXaxis().SetBinLabel(7, "50-100")
# triggerEff.GetXaxis().SetBinLabel(8, "100-1000")
triggerEff.SetLineWidth(2)
triggerEff.SetMaximum(1.0)
triggerEff.SetMinimum(0.2)
c1 = ROOT.TCanvas("c1", "c1", 600, 700)
c1.SetLogx()
triggerEff.Draw("HIST")
latex = ROOT.TLatex()
latex.SetTextAlign(22) 
latex.SetTextSize(0.03)
latex.SetTextFont(42)
latex.SetNDC(True)
latex.DrawLatex(0.5, 0.8, "Iso24 + Mu_17_8 menu")
latex.SetNDC(False)
c1.Update()
c1.SaveAs("TTHH_SL_atleast_2mu_recommended_single_double.png")
