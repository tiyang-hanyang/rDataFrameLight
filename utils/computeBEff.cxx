#include <vector>
#include <string>
#include <fstream>
#include <iostream>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TROOT.h>
#include <TH1D.h>

#include "SampleControl.h"

#include "external/json.hpp"

void getEff(const std::string& era, const std::vector<std::string>& samplePath, const std::string& sampleType)
{
    ROOT::RDataFrame rdf("Events", samplePath);
    ROOT::RDF::RNode rnd(rdf);

    // eta range 0-2.5, not binning, pt range recommended
    Double_t ptBins[] = {20, 30, 50, 70, 100, 140, 200, 300, 600, 1000};
    // flavour has udsg 0, c 4, b 5
    Double_t flaBins[] = {-0.5, 0.5, 4.5, 5.5};

    auto den = rnd.Filter("(JVMweight>0)").Define("GoodJet_pt", "Jet_pt_JEC[GoodJetCond]").Define("GoodJet_hadronFlavour", "Jet_hadronFlavour[GoodJetCond]").Define("totalWeight","genWeight*PUWeight").Histo2D(ROOT::RDF::TH2DModel("den",";p_{T} [GeV];fla", 9, ptBins, 3, flaBins), "GoodJet_pt", "GoodJet_hadronFlavour", "totalWeight");
    std::string denName = sampleType+"s_den";
    TH2D *denClone = (TH2D *)den.GetPtr()->Clone(denName.c_str());

    // add different efficiency measurement condition
    auto rndBtag = rnd.Filter("(JVMweight>0)");
    if (era.find("2024")!=std::string::npos)
        rndBtag = rndBtag.Define("BJet_pt", "Jet_pt_JEC[GoodJetCond && (Jet_btagUParTAK4B > 0.1272)]").Define("BJet_hadronFlavour", "Jet_hadronFlavour[GoodJetCond && (Jet_btagUParTAK4B > 0.1272)]").Define("totalWeight","genWeight*PUWeight");
    else if (era == "Run3Summer23NanoAODv12") 
        rndBtag = rndBtag.Define("BJet_pt", "Jet_pt_JEC[GoodJetCond && (Jet_btagPNetB > 0.1917)]").Define("BJet_hadronFlavour", "Jet_hadronFlavour[GoodJetCond && (Jet_btagPNetB > 0.1917)]").Define("totalWeight","genWeight*PUWeight");
    else if (era == "Run3Summer23BPixNanoAODv12")
    {
        rndBtag = rndBtag.Define("BJet_pt", "Jet_pt_JEC[GoodJetCond && (Jet_btagPNetB > 0.1919)]").Define("BJet_hadronFlavour", "Jet_hadronFlavour[GoodJetCond && (Jet_btagPNetB > 0.1919)]").Define("totalWeight","genWeight*PUWeight");
    }
    else if (era == "Run3Summer22NanoAODv12")
        rndBtag = rndBtag.Define("BJet_pt", "Jet_pt_JEC[GoodJetCond && (Jet_btagPNetB > 0.245)]").Define("BJet_hadronFlavour", "Jet_hadronFlavour[GoodJetCond && (Jet_btagPNetB > 0.245)]").Define("totalWeight","genWeight*PUWeight");
    else if (era == "Run3Summer22EENanoAODv12")
        rndBtag = rndBtag.Define("BJet_pt", "Jet_pt_JEC[GoodJetCond && (Jet_btagPNetB > 0.2605)]").Define("BJet_hadronFlavour", "Jet_hadronFlavour[GoodJetCond && (Jet_btagPNetB > 0.2605)]").Define("totalWeight","genWeight*PUWeight");
    else
    {
        std::cout << "un-recognizaed campaign, please check " << era << std::endl;
        exit(1);
    }

    auto num = rndBtag.Histo2D(ROOT::RDF::TH2DModel("num",";p_{T} [GeV];fla", 9, ptBins, 3, flaBins), "BJet_pt", "BJet_hadronFlavour", "totalWeight");

    std::string numName = sampleType+"_num";
    TH2D *numClone = (TH2D *)num.GetPtr()->Clone(numName.c_str());

    std::string effName = sampleType+"_eff";
    TH2D* eff = (TH2D*)numClone->Clone(effName.c_str());
    eff->Divide(denClone);

    // test output first
    std::vector<double> usdgBinContent; 
    std::vector<double> cBinContent; 
    std::vector<double> bBinContent; 
    // need to first get rid of negative weights
    for (int binIndex=1; binIndex<10; binIndex++)
    {
        auto usdgVal = eff->GetBinContent(binIndex, 1);
        if (usdgVal>1) usdgVal=1;
        if (usdgVal<0) usdgVal=0;
        usdgBinContent.push_back(usdgVal);
        auto cVal = eff->GetBinContent(binIndex, 2);
        if (cVal>1) cVal=1;
        if (cVal<0) cVal=0;
        cBinContent.push_back(cVal);
        auto bVal = eff->GetBinContent(binIndex, 3);
        if (bVal>1) bVal=1;
        if (bVal<0) bVal=0;
        bBinContent.push_back(bVal);
    }

    // write into json
    nlohmann::json saveEff;
    saveEff["schema_version"] = 2;
    saveEff["description"] = "Efficiency for UParTAK4 in 2024_Summer24.";

    // main correction block 
    nlohmann::json correction;
    correction["name"] = "UParTAK4_eff_values";
    correction["description"] = "UParTAK4 efficiency for 2024_Summer24 for b-tagging.";
    correction["version"] = 1;
    nlohmann::json inputBlock = nlohmann::json::parse(R"(
        [
            {
                "name": "systematic", 
                "type": "string"
            },
            {
                "name": "working_point",
                "type": "string",
                "description": "M"
            },
            {
                "name": "flavor", 
                "type": "int",
                "description": "hadron flavor definition: 5=b, 4=c, 0=udsg"
            },
            {
                "name": "abseta", 
                "type": "real"
            },
            {
                "name": "pt", 
                "type": "real"
            }
        ]
    )");
    correction["inputs"] = inputBlock;
    nlohmann::json outputBlock = nlohmann::json::parse(R"(
            {
                "name": "weight",
                "type": "real"
            }
    )");
    correction["output"] = outputBlock;
    // correction block data block
    // correction content usdg
    nlohmann::json udsgPtPartition;
    udsgPtPartition["edges"] = ptBins;
    udsgPtPartition["content"] = usdgBinContent;
    udsgPtPartition["nodetype"] = "binning";
    udsgPtPartition["input"] = "pt";
    udsgPtPartition["flow"] = "clamp";
    nlohmann::json udsgEtaPartition;
    udsgEtaPartition["nodetype"] = "binning";
    udsgEtaPartition["input"] = "abseta";
    udsgEtaPartition["edges"] = nlohmann::json::array({0.0, 2.5});
    udsgEtaPartition["content"] = nlohmann::json::array({udsgPtPartition});
    udsgEtaPartition["flow"] = "error";
    nlohmann::json udsgMWP;
    udsgMWP["key"] = 0;
    udsgMWP["value"] = udsgEtaPartition;
    // correction content c
    nlohmann::json cPtPartition;
    cPtPartition["edges"] = ptBins;
    cPtPartition["content"] = cBinContent;
    cPtPartition["nodetype"] = "binning";
    cPtPartition["input"] = "pt";
    cPtPartition["flow"] = "clamp";
    nlohmann::json cEtaPartition;
    cEtaPartition["nodetype"] = "binning";
    cEtaPartition["input"] = "abseta";
    cEtaPartition["edges"] = nlohmann::json::array({0.0, 2.5});
    cEtaPartition["content"] = nlohmann::json::array({cPtPartition});
    cEtaPartition["flow"] = "error";
    nlohmann::json cMWP;
    cMWP["key"] = 4;
    cMWP["value"] = cEtaPartition;
    // correction content b
    nlohmann::json bPtPartition;
    bPtPartition["edges"] = ptBins;
    bPtPartition["content"] = bBinContent;
    bPtPartition["nodetype"] = "binning";
    bPtPartition["input"] = "pt";
    bPtPartition["flow"] = "clamp";
    nlohmann::json bEtaPartition;
    bEtaPartition["nodetype"] = "binning";
    bEtaPartition["input"] = "abseta";
    bEtaPartition["edges"] = nlohmann::json::array({0.0, 2.5});
    bEtaPartition["content"] = nlohmann::json::array({bPtPartition});
    bEtaPartition["flow"] = "error";
    nlohmann::json bMWP;
    bMWP["key"] = 5;
    bMWP["value"] = bEtaPartition;
    // merging into all the list inside medium WP content
    nlohmann::json mediumContentBlock;
    mediumContentBlock["nodetype"] = "category";
    mediumContentBlock["input"] = "flavor";
    mediumContentBlock["content"] = nlohmann::json::array({udsgMWP, cMWP, bMWP});
    nlohmann::json mediumBlock;
    mediumBlock["key"] = "M";
    mediumBlock["value"] = mediumContentBlock;
    nlohmann::json centralContentBlock;
    centralContentBlock["nodetype"] = "category";
    centralContentBlock["input"] = "working_point";
    centralContentBlock["content"] = nlohmann::json::array({mediumBlock});
    nlohmann::json dataContentBlock;
    dataContentBlock["key"] = "central";
    dataContentBlock["value"] = centralContentBlock;
    nlohmann::json dataBlock;
    dataBlock["nodetype"] = "category";
    dataBlock["input"] = "systematic";
    dataBlock["content"] = nlohmann::json::array({dataContentBlock});
    correction["data"] = dataBlock;
    // add the correction block
    saveEff["corrections"] = nlohmann::json::array({correction});

    // output configs
    std::ofstream o(sampleType+"_btag_eff.json");
    o << saveEff << std::endl;
}

int main(int argc, char** argv)
{
    std::string era = "Run3Summer22EENanoAODv12";
    if (argc >1)
    {
        era = argv[1];
    }
    std::cout << "processing " << era << std::endl;

    ROOT::EnableImplicitMT();
    std::vector<std::string> allProcesses{
        "ttbarDL",
        "ttbarSL",
        "ttbar4Q",

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

        "TbarWp2L",
        "TWm2L",
        //"TbarWp1L",
        //"TWm1L",
        "TZQB",

        "TTHBB",
        "TTHnonBB",

        "WW2L2Nu",
        "WZ2L2Q",
        "ZZ2L2Q",
        "ZZ2L2Nu",
        "WZ3l",
        "ZZ4l",

        "WWW",
        "WWZ",
        "WZZ",
        "ZZZ",

        // // for signal
        // "TTHH_DL_2B2W_batch1",
        // "TTHH_DL_2B2W_batch2",
        // "TTHH_DL_2B2W_batch3",
        // "TTHH_SL_2B2W_batch1",
        // "TTHH_SL_2B2W_batch2",

        "TTBB_DL",
        "TTBB_SL",
        "TTBB_4Q",

        "TTZ_low",
        "TTZ_high",
    };

    if (era == "RunIII2024Summer24NanoAODv15")
    {
        allProcesses.push_back("DY2Mu");
        allProcesses.push_back("DY2Mu_low");

        allProcesses.push_back("WJet_1J");
        allProcesses.push_back("WJet_2J");
        allProcesses.push_back("WJet_3J");
        allProcesses.push_back("WJet_4J");
    }
    else
    {
        allProcesses.push_back("DY2L");
        allProcesses.push_back("DY2L_low");

        allProcesses.push_back("WJets");
        allProcesses.push_back("TTW");
        allProcesses.push_back("TTTT");
    }

    // SampleControl samples("json/samples/RunIII2024Summer24NanoAODv15_passJVM.json");
    std::map<std::string, std::string> jsonPath = {
        {"RunIII2024Summer24NanoAODv15", "/home/tiyang/public/rDataFrameLight_update/source/json/samples/FourJet_NanoAOD/RunIII2024Summer24NanoAODv15_fourJet_temp.json"},
        {"Run3Summer23NanoAODv12", "/home/tiyang/public/rDataFrameLight_update/source/json/samples/FourJet_NanoAOD/Run3Summer23NanoAODv12_fourJet_temp.json"},
        {"Run3Summer23BPixNanoAODv12", "/home/tiyang/public/rDataFrameLight_update/source/json/samples/FourJet_NanoAOD/Run3Summer23BPixNanoAODv12_fourJet_temp.json"},
        {"Run3Summer22NanoAODv12", "/home/tiyang/public/rDataFrameLight_update/source/json/samples/FourJet_NanoAOD/Run3Summer22NanoAODv12_fourJet_temp.json"},
        {"Run3Summer22EENanoAODv12", "/home/tiyang/public/rDataFrameLight_update/source/json/samples/FourJet_NanoAOD/Run3Summer22EENanoAODv12_fourJet_temp.json"},
    };

    SampleControl samples(jsonPath.at(era));

    for (const auto& process:  allProcesses) {
        std::cout << "process: " << process << std::endl;
        getEff(era, samples.getFiles(process), process);
    }

    return 0;
}
