#include "Utility.h"

#include "TFile.h"
#include "TH1.h"
#include "TAxis.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

std::vector<std::string> build_hist_var_names(rdfWS_utility::JsonObject jsonConfig)
{
    std::vector<std::string> histVarNames;

    if (jsonConfig.contains("varNames"))
    {
        for (const auto &variable : jsonConfig.at("varNames").get<std::vector<std::string>>())
            histVarNames.push_back(variable);
    }

    if (jsonConfig.contains("varNames_2D"))
    {
        auto rawVarNames2D = jsonConfig.at("varNames_2D").get<nlohmann::json>();
        if (!rawVarNames2D.is_array())
            rdfWS_utility::messageERROR("mergeCollectedHists", "varNames_2D must be a list of [xVar, yVar] pairs.");

        for (const auto &entry : rawVarNames2D)
        {
            if (!entry.is_array() || entry.size() != 2)
                rdfWS_utility::messageERROR("mergeCollectedHists", "Each varNames_2D entry must be [xVar, yVar].");

            histVarNames.push_back(
                entry.at(0).get<std::string>() + "_vs_" + entry.at(1).get<std::string>()
            );
        }
    }

    if (histVarNames.empty())
        rdfWS_utility::messageERROR("mergeCollectedHists", "At least one variable must be provided in varNames or varNames_2D.");

    return histVarNames;
}

std::vector<std::pair<std::string, std::vector<std::string>>> get_merge_targets(
    rdfWS_utility::JsonObject jsonConfig,
    const std::set<std::string> &selectedChannels)
{
    std::vector<std::pair<std::string, std::vector<std::string>>> mergeTargets;
    std::vector<std::string> datasets = jsonConfig.at("datasets");

    for (const auto &selected : selectedChannels)
    {
        if (std::find(datasets.begin(), datasets.end(), selected) == datasets.end())
            rdfWS_utility::messageERROR("mergeCollectedHists", "Selected channel " + selected + " is not a top-level dataset.");
    }

    nlohmann::json needMerge = nlohmann::json::object();
    if (jsonConfig.contains("needMerge"))
        needMerge = jsonConfig.at("needMerge").get<nlohmann::json>();

    if (!needMerge.is_object())
        rdfWS_utility::messageERROR("mergeCollectedHists", "needMerge must be a map from merged channel to component channel list.");

    for (const auto &dataset : datasets)
    {
        if (!selectedChannels.empty() && selectedChannels.find(dataset) == selectedChannels.end())
            continue;
        if (!needMerge.contains(dataset))
            continue;
        if (!needMerge.at(dataset).is_array())
            rdfWS_utility::messageERROR("mergeCollectedHists", "needMerge entry " + dataset + " must be a list.");

        mergeTargets.push_back({
            dataset,
            needMerge.at(dataset).get<std::vector<std::string>>()
        });
    }

    return mergeTargets;
}

bool axis_compatible(const TAxis *lhs, const TAxis *rhs)
{
    if (!lhs || !rhs)
        return lhs == rhs;
    if (lhs->GetNbins() != rhs->GetNbins())
        return false;

    const TArrayD *lhsBins = lhs->GetXbins();
    const TArrayD *rhsBins = rhs->GetXbins();
    if (lhsBins->GetSize() != rhsBins->GetSize())
        return false;

    if (lhsBins->GetSize() > 0)
    {
        for (int i = 0; i < lhsBins->GetSize(); ++i)
        {
            if (std::fabs(lhsBins->At(i) - rhsBins->At(i)) > 1e-9)
                return false;
        }
        return true;
    }

    return std::fabs(lhs->GetXmin() - rhs->GetXmin()) <= 1e-9
        && std::fabs(lhs->GetXmax() - rhs->GetXmax()) <= 1e-9;
}

bool hist_compatible(const TH1 *lhs, const TH1 *rhs)
{
    if (!lhs || !rhs)
        return false;
    if (std::string(lhs->ClassName()) != std::string(rhs->ClassName()))
        return false;
    if (lhs->GetDimension() != rhs->GetDimension())
        return false;
    if (!axis_compatible(lhs->GetXaxis(), rhs->GetXaxis()))
        return false;
    if (lhs->GetDimension() >= 2 && !axis_compatible(lhs->GetYaxis(), rhs->GetYaxis()))
        return false;
    if (lhs->GetDimension() >= 3 && !axis_compatible(lhs->GetZaxis(), rhs->GetZaxis()))
        return false;
    return true;
}

std::set<std::string> list_hist_suffixes(
    const std::string &componentDir,
    const std::string &component,
    const std::vector<std::string> &histVarNames)
{
    std::set<std::string> suffixes;
    if (!std::filesystem::exists(componentDir))
        return suffixes;

    const std::string prefix = component + "_";
    for (const auto &entry : std::filesystem::directory_iterator(componentDir))
    {
        if (!entry.is_regular_file() || entry.path().extension() != ".root")
            continue;

        const std::string fileName = entry.path().filename().string();
        if (fileName.rfind(prefix, 0) != 0)
            continue;

        const std::string suffix = fileName.substr(prefix.size(), fileName.size() - prefix.size() - 5);
        for (const auto &histVarName : histVarNames)
        {
            if (suffix.rfind(histVarName + "_", 0) == 0)
            {
                suffixes.insert(suffix);
                break;
            }
        }
    }

    return suffixes;
}

std::unique_ptr<TH1> load_hist_clone(const std::string &filePath)
{
    TFile inputFile(filePath.c_str(), "READ");
    if (inputFile.IsZombie())
        rdfWS_utility::messageERROR("mergeCollectedHists", "Failed to open input ROOT file: " + filePath);

    const std::string histName = std::filesystem::path(filePath).stem().string();
    TH1 *inputHist = dynamic_cast<TH1 *>(inputFile.Get(histName.c_str()));
    if (!inputHist)
        rdfWS_utility::messageERROR("mergeCollectedHists", "Failed to get histogram " + histName + " from " + filePath);

    TH1 *histClone = static_cast<TH1 *>(inputHist->Clone(histName.c_str()));
    histClone->SetDirectory(0);
    return std::unique_ptr<TH1>(histClone);
}

void merge_hist_files(
    const std::string &mergedName,
    const std::string &suffix,
    const std::vector<std::string> &inputFiles,
    const std::string &outputDir)
{
    if (inputFiles.empty())
        return;

    auto mergedHist = load_hist_clone(inputFiles[0]);
    for (size_t i = 1; i < inputFiles.size(); ++i)
    {
        auto inputHist = load_hist_clone(inputFiles[i]);
        if (!hist_compatible(mergedHist.get(), inputHist.get()))
        {
            rdfWS_utility::messageERROR(
                "mergeCollectedHists",
                "Histogram shape mismatch while merging " + mergedName + "_" + suffix
                + ": " + inputFiles[0] + " vs " + inputFiles[i]
            );
        }
        mergedHist->Add(inputHist.get());
    }

    rdfWS_utility::creatingFolder("mergeCollectedHists", outputDir);
    const std::string outputHistName = mergedName + "_" + suffix;
    const std::string outputFilePath = outputDir + "/" + outputHistName + ".root";

    TFile outputFile(outputFilePath.c_str(), "RECREATE");
    if (outputFile.IsZombie())
        rdfWS_utility::messageERROR("mergeCollectedHists", "Failed to create output ROOT file: " + outputFilePath);

    mergedHist->SetName(outputHistName.c_str());
    mergedHist->Write("");
}

void merge_channel_hists(
    const std::string &baseDir,
    const std::string &mergedName,
    const std::vector<std::string> &components,
    const std::vector<std::string> &histVarNames)
{
    const std::string componentDir = baseDir + "/" + mergedName;
    std::set<std::string> allSuffixes;
    for (const auto &component : components)
    {
        auto componentSuffixes = list_hist_suffixes(componentDir, component, histVarNames);
        allSuffixes.insert(componentSuffixes.begin(), componentSuffixes.end());
    }

    if (allSuffixes.empty())
    {
        rdfWS_utility::messageWARN("mergeCollectedHists", "Merged channel " + mergedName + " has no component histograms to merge in " + componentDir);
        return;
    }

    for (const auto &suffix : allSuffixes)
    {
        std::vector<std::string> inputFiles;
        for (const auto &component : components)
        {
            const std::string inputFile = componentDir + "/" + component + "_" + suffix + ".root";
            if (!std::filesystem::exists(inputFile))
            {
                rdfWS_utility::messageWARN("mergeCollectedHists", "Missing component histogram, skip this input: " + inputFile);
                continue;
            }
            inputFiles.push_back(inputFile);
        }

        if (inputFiles.empty())
            continue;

        merge_hist_files(mergedName, suffix, inputFiles, baseDir);
    }
}

int main(int argc, char *argv[])
{
    if (argc < 2)
        rdfWS_utility::messageERROR("mergeCollectedHists", "No collect job json provided!");

    const std::string jsonPath = argv[1];
    std::set<std::string> selectedChannels;
    for (int i = 2; i < argc; ++i)
    {
        const std::string arg = argv[i];
        if (arg == "--channels")
            continue;
        selectedChannels.insert(arg);
    }

    rdfWS_utility::JsonObject jsonConfig(rdfWS_utility::readJson("mergeCollectedHists", jsonPath), "JO Config");
    if (jsonConfig.contains("jobType"))
    {
        const std::string jobType = jsonConfig.at("jobType").get<std::string>();
        if (jobType != "collectingHists")
            rdfWS_utility::messageERROR("mergeCollectedHists", "The jobType of your config is not collectingHists! Please check again. Running ceases...");
    }

    const std::vector<std::string> histVarNames = build_hist_var_names(jsonConfig);
    const auto mergeTargets = get_merge_targets(jsonConfig, selectedChannels);
    if (mergeTargets.empty())
    {
        rdfWS_utility::messageWARN("mergeCollectedHists", "No merged channels requested by datasets/needMerge/--channels. Nothing to do.");
        return 0;
    }

    const std::string era = jsonConfig.at("era");
    const std::string baseDir = jsonConfig.at("outDir").get<std::string>() + "_" + era;

    for (const auto &[mergedName, components] : mergeTargets)
    {
        rdfWS_utility::messageINFO("mergeCollectedHists", "Merging channel " + mergedName);
        merge_channel_hists(baseDir, mergedName, components, histVarNames);
    }

    return 0;
}
