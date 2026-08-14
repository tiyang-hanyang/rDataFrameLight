#include "Utility.h"
#include "SkimControl.h"
#include "CutControl.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"

#include <filesystem>
#include <unordered_map>
#include <algorithm>
#include <cctype>

#include "TFile.h"
#include "TH1D.h"

namespace
{
std::function<ROOT::RDF::RNode(ROOT::RDF::RNode)> makeGoldenJsonLambda(const std::string &year)
{
    std::string goldenJsonPath = "json/goldenJson/" + year + "/goldenJson.json";
    auto goldenJson = rdfWS_utility::readJson("skimSamples_oneFile", goldenJsonPath);
    std::map<std::string, std::vector<std::pair<unsigned int, unsigned int>>> goldenJsonRaw = goldenJson;
    std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, unsigned int>>> goldenJsonList;
    goldenJsonList.reserve(goldenJsonRaw.size());

    auto is_digits = [](const std::string &s) {
        return !s.empty() && std::all_of(s.begin(), s.end(), [](unsigned char c) { return std::isdigit(c) != 0; });
    };

    for (const auto &[k, v] : goldenJsonRaw)
    {
        if (!is_digits(k))
        {
            rdfWS_utility::messageWARN("skimSamples_oneFile", "Invalid run key: " + k + ", skip.");
            continue;
        }
        unsigned int run = std::stoul(k);
        goldenJsonList.emplace(run, v);
    }

    return [goldenJsonList](ROOT::RDF::RNode origData) {
        return origData.Filter(
            [goldenJsonList](unsigned int run, unsigned int luminosityBlock) {
                auto it = goldenJsonList.find(run);
                if (it == goldenJsonList.end())
                    return false;
                int isGood = 0;
                for (const auto &blkPair : it->second)
                {
                    if (luminosityBlock <= blkPair.second && luminosityBlock >= blkPair.first)
                    {
                        isGood = 1;
                        break;
                    }
                }
                return isGood != 0;
            },
            {"run", "luminosityBlock"});
    };
}

std::vector<std::string> buildBranchArray(ROOT::RDF::RNode rndDS, const std::vector<std::string> &branchList)
{
    std::vector<std::string> branchArray;
    std::vector<std::string> originalBRs = rndDS.GetColumnNames();
    for (const auto &brName : branchList)
    {
        if (std::find(originalBRs.begin(), originalBRs.end(), brName) == originalBRs.end())
            continue;
        branchArray.push_back(brName);
    }
    return branchArray;
}
} // namespace

int main(int argc, char *argv[])
{
    std::string skimConfig;
    std::string inputFile;
    std::string outputFile;

    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        if (arg == "--skimConfig")
        {
            if (i + 1 >= argc)
            {
                rdfWS_utility::messageERROR("skimSamples_oneFile", "Missing value for --skimConfig");
                return 1;
            }
            skimConfig = argv[++i];
        }
        else if (arg == "--input")
        {
            if (i + 1 >= argc)
            {
                rdfWS_utility::messageERROR("skimSamples_oneFile", "Missing value for --input");
                return 1;
            }
            inputFile = argv[++i];
        }
        else if (arg == "--outFile")
        {
            if (i + 1 >= argc)
            {
                rdfWS_utility::messageERROR("skimSamples_oneFile", "Missing value for --outFile");
                return 1;
            }
            outputFile = argv[++i];
        }
        else
        {
            rdfWS_utility::messageWARN("skimSamples_oneFile",
                                       "Unknown argument: " + arg + ". Expected: --skimConfig --input --outFile");
            return 1;
        }
    }

    if (skimConfig.empty() || inputFile.empty() || outputFile.empty())
    {
        rdfWS_utility::messageERROR(
            "skimSamples_oneFile",
            "Input template wrong. Example: \"skimSamples_oneFile --skimConfig <skim json> --input <input file> --outFile <output file>\"");
        return 1;
    }

    rdfWS_utility::messageINFO("skimSamples_oneFile", "Begin skimming single file");

    auto jsonData = rdfWS_utility::readJson("skimSamples_oneFile", skimConfig);
    rdfWS_utility::JsonObject configFile(jsonData, "Skim JO config");

    if (configFile.contains("jobType"))
    {
        std::string jobType = configFile.at("jobType").get<std::string>();
        if (jobType != "skim")
        {
            rdfWS_utility::messageERROR("skimSamples_oneFile", "The jobType of your config is not skimming!");
            return 1;
        }
    }

    std::string year = configFile.at("year").get<std::string>();
    int isData = configFile.at("isData").get<int>();
    int preliminary = 1;
    if (configFile.contains("preliminary"))
    {
        if (configFile.at("preliminary").get<int>() == 0)
            preliminary = 0;
    }

    std::vector<std::string> cutConfigList = configFile.at("cutConfig").get<std::vector<std::string>>();
    CutControl skimCut;
    if (!cutConfigList.empty())
    {
        skimCut = CutControl(cutConfigList[0]);
        for (size_t i = 1; i < cutConfigList.size(); ++i)
        {
            skimCut = skimCut + CutControl(cutConfigList[i]);
        }
    }

    std::string branchPath = configFile.at("branchConfig").get<std::string>();
    rdfWS_utility::JsonObject branchJson(rdfWS_utility::readJson("skimSamples_oneFile", branchPath), "Branch Config");
    std::vector<std::string> branchList = branchJson.get<std::vector<std::string>>();

    const std::filesystem::path outPath(outputFile);
    if (!outPath.parent_path().empty() && !std::filesystem::exists(outPath.parent_path()))
    {
        std::filesystem::create_directories(outPath.parent_path());
    }
    if (std::filesystem::exists(outputFile))
    {
        TFile fcheck(outputFile.c_str(), "READ");
        if (fcheck.IsZombie() || !fcheck.Get("genWeightSum"))
        {
            const bool isZombie = fcheck.IsZombie();
            fcheck.Close();
            auto baseBackup = outputFile + ".bad";
            std::string backupPath = baseBackup;
            int suffix = 0;
            while (std::filesystem::exists(backupPath))
            {
                ++suffix;
                backupPath = baseBackup + "." + std::to_string(suffix);
            }
            rdfWS_utility::messageWARN("skimSamples_oneFile",
                                       std::string("Output ") + (isZombie ? "is zombie" : "missing genWeightSum") +
                                           ", backing up to: " + backupPath);
            std::filesystem::rename(outputFile, backupPath);
        }
        else
        {
            rdfWS_utility::messageINFO("skimSamples_oneFile", outputFile + " already exists and is valid, skip");
            return 0;
        }
    }

    ROOT::RDataFrame rdfDS("Events", inputFile.c_str());
    ROOT::RDF::RNode rndDS(rdfDS);

    ROOT::RDF::RResultPtr<float> sumw;
    if (!isData)
    {
        sumw = rndDS.Sum<float>("genWeight");
    }
    else
    {
        sumw = rndDS.Define("dummyWeight", "1.0f").Sum<float>("dummyWeight");
    }

    if (preliminary && isData)
    {
        auto goldenJsonLambda = makeGoldenJsonLambda(year);
        rndDS = goldenJsonLambda(rndDS);
    }

    if (!cutConfigList.empty())
    {
        rndDS = skimCut.applyCut(rndDS);
    }

    auto branchArray = buildBranchArray(rndDS, branchList);

    ROOT::RDF::RSnapshotOptions SSoption;
    SSoption.fCompressionLevel = 4;
    SSoption.fLazy = true;

    auto snap = rndDS.Snapshot("Events", outputFile, branchArray, SSoption);
    std::vector<ROOT::RDF::RResultHandle> handles;
    handles.reserve(2);
    handles.emplace_back(ROOT::RDF::RResultHandle(snap));
    handles.emplace_back(ROOT::RDF::RResultHandle(sumw));
    ROOT::RDF::RunGraphs(handles);

    {
        TFile fRC(outputFile.c_str(), "UPDATE");
        if (!fRC.IsOpen())
        {
            rdfWS_utility::messageWARN("skimSamples_oneFile", "Failed to open output file for genWeightSum: " + outputFile);
            return 1;
        }
        auto h_sumw = new TH1D("genWeightSum", "sum of genWeight (this file)", 1, 0.0, 1.0);
        h_sumw->SetBinContent(1, static_cast<double>(sumw.GetValue()));
        h_sumw->GetYaxis()->SetTitle("sum(genWeight)");
        h_sumw->Write("");
        fRC.Close();
    }

    return 0;
}
