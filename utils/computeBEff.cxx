#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include <ROOT/RDataFrame.hxx>

#include "SampleControl.h"

#include "external/json.hpp"

double get_2024_wp_threshold(const std::string& workingPoint)
{
    if (workingPoint == "L") return 0.0246;
    if (workingPoint == "M") return 0.1272;
    if (workingPoint == "T") return 0.4648;
    throw std::runtime_error("Unknown 2024 working point: " + workingPoint);
}

std::string build_tagged_jet_expr(const std::string& jetMaskExpr, const std::string& discriminatorCut)
{
    return jetMaskExpr + " && (" + discriminatorCut + ")";
}

void print_usage(const char* programName)
{
    std::cout << "Usage: " << programName
              << " ERA SAMPLE_JSON [--channel CHANNEL ...] [--good-jet-expr SPEC] [--muon-pt-baseline-expr EXPR]"
              << std::endl;
    std::cout << "Arguments:" << std::endl;
    std::cout << "  ERA          Required campaign / era name." << std::endl;
    std::cout << "  SAMPLE_JSON  Required sample json path." << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  --channel CHANNEL ..." << std::endl;
    std::cout << "      Optional. After --channel, all following plain arguments are treated" << std::endl;
    std::cout << "      as channel names until the next --option appears. Omit this option" << std::endl;
    std::cout << "      to process all channels from SAMPLE_JSON." << std::endl;
    std::cout << "  --good-jet-expr SPEC" << std::endl;
    std::cout << "      Optional and repeatable. Default is GoodJetCond." << std::endl;
    std::cout << "      SPEC may be an existing branch name or a full RDataFrame expression." << std::endl;
    std::cout << "  --muon-pt-baseline-expr EXPR" << std::endl;
    std::cout << "      Optional. Default is the built-in dimuon baseline using Muon_pt_Rcorr," << std::endl;
    std::cout << "      leadingMuonIdx, and subleadingMuonIdx. Provide a custom expression to" << std::endl;
    std::cout << "      override that baseline, e.g. for one-muon samples." << std::endl;
    std::cout << "Examples:" << std::endl;
    std::cout << "  " << programName << " RunIII2024Summer24NanoAODv15 /path/to/sample.json" << std::endl;
    std::cout << "  " << programName << " RunIII2024Summer24NanoAODv15 /path/to/sample.json --channel ttbarDL TTBB_DL" << std::endl;
    std::cout << "  " << programName << " RunIII2024Summer24NanoAODv15 /path/to/sample.json --good-jet-expr '(Jet_pt_JEC > 30.0) && (Jet_rawFactor < 0.9) && Jet_passJetIdTight'" << std::endl;
    std::cout << "  " << programName << " RunIII2024Summer24NanoAODv15 /path/to/sample.json --muon-pt-baseline-expr '(leadingMuonIdx < Muon_pt_Rcorr.size()) && (Muon_pt_Rcorr[leadingMuonIdx] > 20.0)'" << std::endl;
}

void getEff(
    const std::string& era,
    const std::vector<std::string>& samplePath,
    const std::string& sampleType,
    const std::string& goodJetSpec,
    const std::string& muonPtBaselineExpr
)
{
    std::vector<std::string> treeExistSamples;
    for (const auto& singleSample : samplePath)
    {
        TFile* fTemp = new TFile(singleSample.c_str(), "read");
        if (!fTemp || fTemp->IsZombie())
        {
            std::cout << "file not exist: " << singleSample << std::endl;
            continue;
        }
        if (!fTemp->GetListOfKeys()->FindObject("Events"))
        {
            std::cout << "no Events in " << singleSample << std::endl;
            fTemp->Close();
            continue;
        }
        treeExistSamples.push_back(singleSample);
        fTemp->Close();
    }
    if (treeExistSamples.size() == 0)
    {
        std::cout << "No available files for " << sampleType << std::endl;
        return;
    }

    ROOT::RDataFrame rdf("Events", treeExistSamples);
    ROOT::RDF::RNode rnd(rdf);

    const std::string defaultMuonPtBaseline =
        "(leadingMuonIdx < Muon_pt_Rcorr.size()) && "
        "(subleadingMuonIdx < Muon_pt_Rcorr.size()) && "
        "(Muon_pt_Rcorr[leadingMuonIdx] > 20.0) && "
        "(Muon_pt_Rcorr[subleadingMuonIdx] > 15.0)";
    const bool useDefaultMuonPtBaseline = muonPtBaselineExpr.empty();
    const std::string muonPtBaseline = useDefaultMuonPtBaseline ? defaultMuonPtBaseline : muonPtBaselineExpr;
    if (useDefaultMuonPtBaseline)
    {
        if (!rnd.HasColumn("Muon_pt_Rcorr") || !rnd.HasColumn("leadingMuonIdx") || !rnd.HasColumn("subleadingMuonIdx"))
        {
            throw std::runtime_error(
                "computeBEff requires Muon_pt_Rcorr, leadingMuonIdx, and subleadingMuonIdx to apply the default muon baseline. "
                "Provide --muon-pt-baseline-expr to override it."
            );
        }
        std::cout << "[computeBEff] using default muon baseline: " << muonPtBaseline << std::endl;
    }
    else
    {
        std::cout << "[computeBEff] using custom muon baseline: " << muonPtBaseline << std::endl;
    }
    rnd = rnd.Filter(muonPtBaseline);

    const std::string goodJetTempName = "GoodJetCond_eff_temp";
    const bool isExistingBranch = rnd.HasColumn(goodJetSpec);
    if (isExistingBranch)
    {
        std::cout << "[computeBEff] good jet spec '" << goodJetSpec << "' matches an existing branch." << std::endl;
        std::cout << "[computeBEff] defining alias: " << goodJetTempName << " = " << goodJetSpec << std::endl;
    }
    else
    {
        std::cout << "[computeBEff] good jet spec '" << goodJetSpec << "' does not match an existing branch." << std::endl;
        std::cout << "[computeBEff] defining expression: " << goodJetTempName << " = " << goodJetSpec << std::endl;
    }
    rnd = rnd.Define(goodJetTempName, goodJetSpec);

    Double_t ptBins[] = {20, 30, 50, 70, 100, 140, 200, 300, 600, 1000};
    Double_t flaBins[] = {-0.5, 0.5, 4.5, 5.5};
    nlohmann::json ptEdges = nlohmann::json::array();
    for (int i = 0; i < 10; ++i)
        ptEdges.push_back(ptBins[i]);
    const std::string goodJetMaskExpr = "(" + goodJetTempName + ")";

    auto den = rnd.Filter("(JVMweight>0)")
        .Define("GoodJet_pt", "Jet_pt_JEC[" + goodJetMaskExpr + "]")
        .Define("GoodJet_hadronFlavour", "Jet_hadronFlavour[" + goodJetMaskExpr + "]")
        .Define("totalWeight","genWeight*PUWeight")
        .Histo2D(ROOT::RDF::TH2DModel("den",";p_{T} [GeV];fla", 9, ptBins, 3, flaBins), "GoodJet_pt", "GoodJet_hadronFlavour", "totalWeight");
    std::string denName = sampleType + "s_den";
    TH2D *denClone = (TH2D *)den.GetPtr()->Clone(denName.c_str());

    std::vector<std::string> workingPoints = {"M"};
    if (era.find("2024") != std::string::npos)
        workingPoints = {"L", "M", "T"};

    nlohmann::json workingPointContent = nlohmann::json::array();
    for (const auto& workingPoint : workingPoints)
    {
        auto rndBtag = rnd.Filter("(JVMweight>0)");
        std::string cutExpr;
        if (era.find("2024") != std::string::npos)
        {
            auto threshold = std::to_string(get_2024_wp_threshold(workingPoint));
            cutExpr = build_tagged_jet_expr(goodJetMaskExpr, "Jet_btagUParTAK4B > " + threshold);
        }
        else if (era == "Run3Summer23NanoAODv12")
            cutExpr = build_tagged_jet_expr(goodJetMaskExpr, "Jet_btagPNetB > 0.1917");
        else if (era.find("Run3Summer23BPixNanoAODv12") != std::string::npos)
            cutExpr = build_tagged_jet_expr(goodJetMaskExpr, "Jet_btagPNetB > 0.1919");
        else if (era == "Run3Summer22NanoAODv12")
            cutExpr = build_tagged_jet_expr(goodJetMaskExpr, "Jet_btagPNetB > 0.245");
        else if (era == "Run3Summer22EENanoAODv12")
            cutExpr = build_tagged_jet_expr(goodJetMaskExpr, "Jet_btagPNetB > 0.2605");
        else
        {
            std::cout << "un-recognized campaign, please check " << era << std::endl;
            exit(1);
        }

        rndBtag = rndBtag.Define("BJet_pt", "Jet_pt_JEC[" + cutExpr + "]")
                         .Define("BJet_hadronFlavour", "Jet_hadronFlavour[" + cutExpr + "]")
                         .Define("totalWeight", "genWeight*PUWeight");

        auto num = rndBtag.Histo2D(ROOT::RDF::TH2DModel("num",";p_{T} [GeV];fla", 9, ptBins, 3, flaBins), "BJet_pt", "BJet_hadronFlavour", "totalWeight");

        std::string numName = sampleType + "_" + workingPoint + "_num";
        TH2D *numClone = (TH2D *)num.GetPtr()->Clone(numName.c_str());

        std::string effName = sampleType + "_" + workingPoint + "_eff";
        TH2D* eff = (TH2D*)numClone->Clone(effName.c_str());
        eff->Divide(denClone);

        std::vector<double> usdgBinContent;
        std::vector<double> cBinContent;
        std::vector<double> bBinContent;
        for (int binIndex = 1; binIndex < 10; ++binIndex)
        {
            auto usdgVal = eff->GetBinContent(binIndex, 1);
            if (usdgVal > 1) usdgVal = 1;
            if (usdgVal < 0) usdgVal = 0;
            usdgBinContent.push_back(usdgVal);

            auto cVal = eff->GetBinContent(binIndex, 2);
            if (cVal > 1) cVal = 1;
            if (cVal < 0) cVal = 0;
            cBinContent.push_back(cVal);

            auto bVal = eff->GetBinContent(binIndex, 3);
            if (bVal > 1) bVal = 1;
            if (bVal < 0) bVal = 0;
            bBinContent.push_back(bVal);
        }

        nlohmann::json udsgPtPartition;
        udsgPtPartition["edges"] = ptEdges;
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

        nlohmann::json cPtPartition;
        cPtPartition["edges"] = ptEdges;
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

        nlohmann::json bPtPartition;
        bPtPartition["edges"] = ptEdges;
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

        nlohmann::json udsgEntry;
        udsgEntry["key"] = 0;
        udsgEntry["value"] = udsgEtaPartition;

        nlohmann::json cEntry;
        cEntry["key"] = 4;
        cEntry["value"] = cEtaPartition;

        nlohmann::json bEntry;
        bEntry["key"] = 5;
        bEntry["value"] = bEtaPartition;

        nlohmann::json flavorCategory;
        flavorCategory["nodetype"] = "category";
        flavorCategory["input"] = "flavor";
        flavorCategory["content"] = nlohmann::json::array({udsgEntry, cEntry, bEntry});

        nlohmann::json workingPointEntry;
        workingPointEntry["key"] = workingPoint;
        workingPointEntry["value"] = flavorCategory;
        workingPointContent.push_back(workingPointEntry);
    }

    nlohmann::json saveEff;
    saveEff["schema_version"] = 2;
    saveEff["description"] = "Efficiency for UParTAK4 in 2024_Summer24.";

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

    nlohmann::json centralContentBlock;
    centralContentBlock["nodetype"] = "category";
    centralContentBlock["input"] = "working_point";
    centralContentBlock["content"] = workingPointContent;

    nlohmann::json dataContentBlock;
    dataContentBlock["key"] = "central";
    dataContentBlock["value"] = centralContentBlock;

    nlohmann::json dataBlock;
    dataBlock["nodetype"] = "category";
    dataBlock["input"] = "systematic";
    dataBlock["content"] = nlohmann::json::array({dataContentBlock});
    correction["data"] = dataBlock;

    saveEff["corrections"] = nlohmann::json::array({correction});

    std::ofstream o(sampleType + "_btag_eff.json");
    o << saveEff << std::endl;
}

int main(int argc, char** argv)
{
    if ((argc > 1 && (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h")) || argc < 3)
    {
        print_usage(argv[0]);
        return 1;
    }

    const std::string era = argv[1];
    const std::string jsonPath = argv[2];
    std::cout << "processing era " << era << std::endl;
    std::cout << "using sample json " << jsonPath << std::endl;

    ROOT::EnableImplicitMT();
    SampleControl samples(jsonPath);
    const auto allAvailableChannels = samples.getAllChannels();

    std::vector<std::string> channels;
    std::vector<std::string> goodJetSpecs = {"GoodJetCond"};
    std::string muonPtBaselineExpr;
    for (int i = 3; i < argc; ++i)
    {
        const std::string arg = argv[i];
        if (arg == "--channel")
        {
            if (i + 1 >= argc || std::string(argv[i + 1]).rfind("--", 0) == 0)
                throw std::runtime_error("--channel requires a value.");
            while (i + 1 < argc)
            {
                const std::string nextArg = argv[i + 1];
                if (nextArg.rfind("--", 0) == 0)
                    break;
                channels.push_back(nextArg);
                ++i;
            }
        }
        else if (arg == "--good-jet-expr")
        {
            if (i + 1 >= argc)
                throw std::runtime_error("--good-jet-expr requires a value.");
            if (goodJetSpecs.size() == 1 && goodJetSpecs[0] == "GoodJetCond")
                goodJetSpecs.clear();
            goodJetSpecs.push_back(argv[++i]);
        }
        else if (arg == "--muon-pt-baseline-expr")
        {
            if (i + 1 >= argc)
                throw std::runtime_error("--muon-pt-baseline-expr requires a value.");
            muonPtBaselineExpr = argv[++i];
        }
        else
        {
            throw std::runtime_error("Unknown option or positional argument: " + arg);
        }
    }

    if (channels.empty())
    {
        channels = allAvailableChannels;
    }

    std::cout << "channels to scan: " << channels.size() << std::endl;
    std::cout << "good jet specs to scan: " << goodJetSpecs.size() << std::endl;
    for (const auto& channel : channels)
    {
        std::cout << "process: " << channel << std::endl;
        for (const auto& goodJetSpec : goodJetSpecs)
        {
            std::cout << "good jet spec: " << goodJetSpec << std::endl;
            getEff(era, samples.getFiles(channel), channel, goodJetSpec, muonPtBaselineExpr);
        }
    }

    return 0;
}
