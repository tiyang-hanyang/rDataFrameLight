#include "PlotControl.h"
#include "HistControl.h"

#include "Utility.h"

#include <iostream>
#include <exception>
#include <string>
#include <map>
#include <vector>
#include <sstream>
#include <memory>
#include <cmath>
#include <utility>
#include <algorithm>
#include <set>

#include "TH1D.h"
#include "TFile.h"

// function to read systmatic variables
// if non provided, will return empty vector and plot nominal-only
std::vector<std::map<std::string, double>> getStackUncert(std::vector<std::string> stackOrder, rdfWS_utility::JsonObject jsonConfig)
{
    if (stackOrder.size() == 0)
        return {};
    std::string systFile = "";
    if (jsonConfig.contains("globalSystConfig"))
        systFile = jsonConfig.at("globalSystConfig").get<std::string>();
    else if (jsonConfig.contains("systConfig"))
        systFile = jsonConfig.at("systConfig").get<std::string>();
    if (systFile == "") return {};

    std::map<std::string, std::vector<double>> systJson = rdfWS_utility::readJson("plotHists", systFile);

    std::map<std::string, double> systUp, systDown;
    for (auto &[key, variation] : systJson)
    {
        systUp.emplace(key, variation[1]);
        systDown.emplace(key, variation[0]);
    }
    return {systUp, systDown};
}

namespace
{
TGraphAsymmErrors *buildStandaloneUncertaintyBand(
    TH1D *nominalHist,
    const std::vector<std::pair<TH1D *, TH1D *>> &shapeUncerts,
    bool includeMCStat)
{
    if (nominalHist == nullptr)
        return nullptr;
    if (!includeMCStat && shapeUncerts.empty())
        return nullptr;

    TGraphAsymmErrors *systBand = new TGraphAsymmErrors(nominalHist);
    for (int i = 1; i <= nominalHist->GetNbinsX(); ++i)
    {
        const double binCenter = nominalHist->GetBinCenter(i);
        const double binWidth = nominalHist->GetBinWidth(i) / 2.0;
        const double nominalVal = nominalHist->GetBinContent(i);

        double errUp2 = 0.0;
        double errDown2 = 0.0;
        if (includeMCStat)
        {
            const double statErr = nominalHist->GetBinError(i);
            errUp2 += statErr * statErr;
            errDown2 += statErr * statErr;
        }

        auto addEnvelope = [&](double upVal, double downVal)
        {
            const double upDelta = upVal - nominalVal;
            const double downDelta = downVal - nominalVal;
            const double upErr = std::max({0.0, upDelta, downDelta});
            const double downErr = std::max({0.0, -upDelta, -downDelta});
            errUp2 += upErr * upErr;
            errDown2 += downErr * downErr;
        };

        for (const auto &[shapeUp, shapeDown] : shapeUncerts)
        {
            if (shapeUp == nullptr || shapeDown == nullptr)
                continue;
            addEnvelope(shapeUp->GetBinContent(i), shapeDown->GetBinContent(i));
        }

        systBand->SetPoint(i - 1, binCenter, nominalVal);
        systBand->SetPointError(i - 1, binWidth, binWidth, std::sqrt(errDown2), std::sqrt(errUp2));
    }

    systBand->SetFillColor(nominalHist->GetLineColor());
    systBand->SetFillStyle(3004);
    systBand->SetLineWidth(0);
    systBand->SetMarkerSize(0);
    return systBand;
}

std::vector<int> getCropBinIndicesFromConfig(
    rdfWS_utility::JsonObject varConfig,
    const std::vector<double> &cropRange)
{
    std::vector<int> selectedBins;
    if (cropRange.size() != 2)
        rdfWS_utility::messageERROR("plotHists", "cropedRange entry must have exactly two values [min, max].");

    const double cropMin = cropRange[0];
    const double cropMax = cropRange[1];

    if (varConfig.contains("binning"))
    {
        const auto binEdges = varConfig.at("binning").get<std::vector<double>>();
        if (binEdges.size() < 2)
            rdfWS_utility::messageERROR("plotHists", "Variable-binning config must contain at least two edges.");
        for (size_t i = 0; i + 1 < binEdges.size(); ++i)
        {
            if (binEdges[i] >= cropMin && binEdges[i + 1] <= cropMax)
                selectedBins.push_back(static_cast<int>(i));
        }
        return selectedBins;
    }

    const int origBins = std::stoi(varConfig.at("nBins").get<std::string>());
    const double origMin = std::stod(varConfig.at("min").get<std::string>());
    const double origMax = std::stod(varConfig.at("max").get<std::string>());
    const double binWidth = (origMax - origMin) / origBins;
    for (int i = 0; i < origBins; ++i)
    {
        const double lowEdge = origMin + i * binWidth;
        const double upEdge = lowEdge + binWidth;
        if (lowEdge >= cropMin && upEdge <= cropMax)
            selectedBins.push_back(i);
    }
    return selectedBins;
}

std::vector<std::pair<double, double>> getBlindRanges(
    rdfWS_utility::JsonObject jsonConfig,
    const std::string &varName)
{
    if (!jsonConfig.contains("blindRange"))
        return {};

    auto blindConfig = jsonConfig.at("blindRange").get<nlohmann::json>();
    if (!blindConfig.is_object() || !blindConfig.contains(varName))
        return {};

    const auto rangeConfig = blindConfig.at(varName);
    std::vector<std::pair<double, double>> blindRanges;
    if (rangeConfig.is_array() && rangeConfig.size() == 2 && rangeConfig[0].is_number() && rangeConfig[1].is_number())
    {
        blindRanges.emplace_back(rangeConfig[0].get<double>(), rangeConfig[1].get<double>());
        return blindRanges;
    }
    if (rangeConfig.is_array())
    {
        for (const auto &item : rangeConfig)
        {
            if (!item.is_array() || item.size() != 2)
                rdfWS_utility::messageERROR("plotHists", "blindRange." + varName + " must be [min, max] or a list of such pairs.");
            blindRanges.emplace_back(item[0].get<double>(), item[1].get<double>());
        }
        return blindRanges;
    }

    rdfWS_utility::messageERROR("plotHists", "blindRange." + varName + " must be [min, max] or a list of such pairs.");
    return {};
}

std::map<std::string, std::vector<std::string>> loadSystAliases()
{
    rdfWS_utility::JsonObject aliasConfig(
        rdfWS_utility::readJson("plotHists", "json/general_config/syst_NP.json"),
        "Systematic Alias Config"
    );
    return aliasConfig.get<std::map<std::string, std::vector<std::string>>>();
}

std::vector<std::string> expandAliasOrName(
    const std::string &name,
    const std::map<std::string, std::vector<std::string>> &aliases)
{
    auto iter = aliases.find(name);
    if (iter != aliases.end())
        return iter->second;
    return {name};
}

float sumLumiValues(
    rdfWS_utility::JsonObject lumiConfig,
    const std::vector<std::string> &eras,
    const std::string &context)
{
    if (eras.empty())
        rdfWS_utility::messageERROR("plotHists", "Empty era list for lumi sum: " + context);
    float total = lumiConfig.at(eras[0]).get<float>();
    for (size_t i = 1; i < eras.size(); ++i)
        total += lumiConfig.at(eras[i]).get<float>();
    return total;
}

std::set<std::string> expandSignalLoadChannels(
    const std::vector<std::string> &signalChannels,
    const std::map<std::string, std::vector<std::string>> &needMerge)
{
    std::set<std::string> expanded;
    for (const auto &channel : signalChannels)
    {
        expanded.insert(channel);
        auto mergeIter = needMerge.find(channel);
        if (mergeIter != needMerge.end())
        {
            for (const auto &component : mergeIter->second)
                expanded.insert(component);
        }
    }
    return expanded;
}

std::string summarizeJsonValue(const nlohmann::json &value, size_t maxLen = 200)
{
    std::string dumped = value.dump();
    if (dumped.size() > maxLen)
        dumped = dumped.substr(0, maxLen) + "...";
    return dumped;
}
}

struct MaxVariationLinePlotConfig
{
    bool enabled = false;
    std::string channel;
};

std::map<std::string, std::vector<std::string>> getWeightSystematics(
    rdfWS_utility::JsonObject jsonConfig,
    const std::string &fieldName = "weightSyst")
{
    if (!jsonConfig.contains(fieldName))
        return {};
    const auto aliases = loadSystAliases();
    std::map<std::string, std::vector<std::string>> expanded;
    auto rawConfig = jsonConfig.at(fieldName).get<nlohmann::json>();
    if (!rawConfig.is_object())
        rdfWS_utility::messageERROR("plotHists", fieldName + " must be a map.");

    for (auto it = rawConfig.begin(); it != rawConfig.end(); ++it)
    {
        const std::string rawName = it.key();
        if (!it.value().is_array())
            rdfWS_utility::messageERROR("plotHists", fieldName + " entry " + rawName + " must be a list.");
        const auto config = it.value().get<std::vector<std::string>>();
        const auto expandedNames = expandAliasOrName(rawName, aliases);
        if (expandedNames.size() > 1 && config.size() != 1)
            rdfWS_utility::messageERROR("plotHists", fieldName + " alias " + rawName + " only supports shorthand form [nominalWeight].");
        for (const auto &name : expandedNames)
        {
            if (config.size() == 1)
            {
                expanded.emplace(name, std::vector<std::string>{
                    config[0],
                    config[0] + "_" + name + "_up",
                    config[0] + "_" + name + "_down"
                });
            }
            else
            {
                expanded.emplace(name, config);
            }
        }
    }
    return expanded;
}

struct OverlaySystConfig
{
    std::vector<std::string> channels;
    bool includeMCStat = true;
};

OverlaySystConfig getOverlaySystConfig(
    rdfWS_utility::JsonObject jsonConfig,
    const std::vector<std::string> &channels,
    const std::vector<std::string> &stackOrder,
    const std::map<std::string, int> &isData,
    const std::map<std::string, int> &isDD)
{
    OverlaySystConfig config;
    if (!jsonConfig.contains("overlaySyst"))
        return config;

    auto overlayConfig = jsonConfig.at("overlaySyst");
    if (overlayConfig.contains("channels"))
        config.channels = overlayConfig.at("channels").get<std::vector<std::string>>();
    if (overlayConfig.contains("includeMCStat"))
        config.includeMCStat = overlayConfig.at("includeMCStat").get<bool>();

    if (config.channels.empty())
    {
        for (const auto &channel : channels)
        {
            auto isDataIter = isData.find(channel);
            if (isDataIter != isData.end() && isDataIter->second == 1)
                continue;
            auto isDDIter = isDD.find(channel);
            if (isDDIter != isDD.end() && isDDIter->second == 1)
                continue;
            if (std::find(stackOrder.begin(), stackOrder.end(), channel) != stackOrder.end())
                continue;
            config.channels.push_back(channel);
        }
    }

    return config;
}

std::vector<std::string> getShiftSystematicNames(
    rdfWS_utility::JsonObject jsonConfig,
    const std::string &shiftFieldName = "shiftSyst",
    const std::string &shiftExprFieldName = "shiftExprSyst")
{
    std::set<std::string> shiftSystSet;
    const auto aliases = loadSystAliases();
    if (jsonConfig.contains(shiftFieldName))
    {
        auto shiftConfig = jsonConfig.at(shiftFieldName).get<nlohmann::json>();
        if (!shiftConfig.is_object())
            rdfWS_utility::messageERROR("plotHists", shiftFieldName + " must be a map.");
        for (auto it = shiftConfig.begin(); it != shiftConfig.end(); ++it)
        {
            for (const auto &name : expandAliasOrName(it.key(), aliases))
                shiftSystSet.insert(name);
        }
    }
    if (jsonConfig.contains(shiftExprFieldName))
    {
        auto shiftExprConfig = jsonConfig.at(shiftExprFieldName).get<nlohmann::json>();
        if (!shiftExprConfig.is_object())
            rdfWS_utility::messageERROR("plotHists", shiftExprFieldName + " must be a map.");
        for (auto it = shiftExprConfig.begin(); it != shiftExprConfig.end(); ++it)
        {
            for (const auto &name : expandAliasOrName(it.key(), aliases))
                shiftSystSet.insert(name);
        }
    }
    std::vector<std::string> shiftSysts(shiftSystSet.begin(), shiftSystSet.end());
    return shiftSysts;
}

std::map<std::string, int> getIsDataMap(
    rdfWS_utility::JsonObject jsonConfig,
    const std::vector<std::string> &datasets)
{
    std::map<std::string, int> flags;
    for (const auto &dataset : datasets)
    {
        flags.emplace(dataset, 0);
    }

    if (!jsonConfig.contains("isData"))
        return flags;

    auto rawConfig = jsonConfig.at("isData").get<nlohmann::json>();
    if (rawConfig.is_array())
    {
        for (const auto &item : rawConfig.get<std::vector<std::string>>())
        {
            flags[item] = 1;
        }
        return flags;
    }

    if (rawConfig.is_object())
    {
        for (auto it = rawConfig.begin(); it != rawConfig.end(); ++it)
        {
            if (!it.value().is_boolean() && !it.value().is_number_integer())
                rdfWS_utility::messageERROR(
                    "plotHists",
                    "isData object entry " + it.key() + " must be an integer or boolean flag."
                );
            flags[it.key()] = it.value().get<int>();
        }
        return flags;
    }

    rdfWS_utility::messageERROR(
        "plotHists",
        "isData must be a list of dataset names or a map of dataset flags, but got type "
        + std::string(rawConfig.type_name()) + " with value " + summarizeJsonValue(rawConfig)
    );

    return flags;
}

std::map<std::string, int> getIsDDMap(
    rdfWS_utility::JsonObject jsonConfig,
    const std::vector<std::string> &datasets)
{
    std::map<std::string, int> flags;
    for (const auto &dataset : datasets)
    {
        flags.emplace(dataset, 0);
    }

    if (!jsonConfig.contains("isDD"))
        return flags;

    auto rawConfig = jsonConfig.at("isDD").get<nlohmann::json>();
    if (rawConfig.is_array())
    {
        for (const auto &item : rawConfig.get<std::vector<std::string>>())
        {
            flags[item] = 1;
        }
        return flags;
    }

    if (rawConfig.is_object())
    {
        for (auto it = rawConfig.begin(); it != rawConfig.end(); ++it)
        {
            if (!it.value().is_boolean() && !it.value().is_number_integer())
                rdfWS_utility::messageERROR(
                    "plotHists",
                    "isDD object entry " + it.key() + " must be an integer or boolean flag."
                );
            flags[it.key()] = it.value().get<int>();
        }
        return flags;
    }

    rdfWS_utility::messageERROR(
        "plotHists",
        "isDD must be a list of dataset names or a map of dataset flags, but got type "
        + std::string(rawConfig.type_name()) + " with value " + summarizeJsonValue(rawConfig)
    );

    return flags;
}

std::map<std::string, double> getDDVariationMap(
    rdfWS_utility::JsonObject jsonConfig,
    const std::vector<std::string> &datasets)
{
    std::map<std::string, double> variations;
    for (const auto &dataset : datasets)
        variations.emplace(dataset, 0.0);

    std::string fieldName = "";
    if (jsonConfig.contains("DDvariation"))
        fieldName = "DDvariation";
    else if (jsonConfig.contains("DDvartion"))
        fieldName = "DDvartion";

    if (fieldName == "")
        return variations;

    auto rawConfig = jsonConfig.at(fieldName).get<nlohmann::json>();
    if (!rawConfig.is_object())
    {
        rdfWS_utility::messageERROR(
            "plotHists",
            fieldName + " must be a map from DD dataset name to relative uncertainty."
        );
    }

    for (auto it = rawConfig.begin(); it != rawConfig.end(); ++it)
    {
        if (!it.value().is_number())
        {
            rdfWS_utility::messageERROR(
                "plotHists",
                fieldName + " entry " + it.key() + " must be a numeric relative uncertainty."
            );
        }
        variations[it.key()] = it.value().get<double>();
    }

    return variations;
}

std::map<std::string, double> getProcessScaleMap(
    rdfWS_utility::JsonObject jsonConfig,
    const std::vector<std::string> &datasets)
{
    std::map<std::string, double> scales;
    for (const auto &dataset : datasets)
        scales.emplace(dataset, 1.0);

    if (!jsonConfig.contains("processScale"))
        return scales;

    auto rawConfig = jsonConfig.at("processScale").get<nlohmann::json>();
    if (!rawConfig.is_object())
    {
        rdfWS_utility::messageERROR(
            "plotHists",
            "processScale must be a map from dataset/component name to multiplicative scale."
        );
    }

    for (auto it = rawConfig.begin(); it != rawConfig.end(); ++it)
    {
        if (!it.value().is_number())
        {
            rdfWS_utility::messageERROR(
                "plotHists",
                "processScale entry " + it.key() + " must be numeric."
            );
        }
        scales[it.key()] = it.value().get<double>();
    }
    return scales;
}

MaxVariationLinePlotConfig getMaxVariationLinePlotConfig(
    rdfWS_utility::JsonObject jsonConfig,
    const std::vector<std::string> &datasets,
    const std::map<std::string, int> &isData,
    const std::map<std::string, int> &isDD)
{
    MaxVariationLinePlotConfig config;
    if (!jsonConfig.contains("plotMode"))
        return config;

    const std::string plotMode = jsonConfig.at("plotMode").get<std::string>();
    if (plotMode != "maxVariationLines")
        return config;

    config.enabled = true;
    if (jsonConfig.contains("maxVariationChannel"))
    {
        config.channel = jsonConfig.at("maxVariationChannel").get<std::string>();
    }
    else
    {
        std::vector<std::string> nonDataChannels;
        for (const auto &dataset : datasets)
        {
            auto isDataIter = isData.find(dataset);
            if (isDataIter != isData.end() && isDataIter->second == 1)
                continue;
            auto isDDIter = isDD.find(dataset);
            if (isDDIter != isDD.end() && isDDIter->second == 1)
                continue;
            nonDataChannels.push_back(dataset);
        }
        if (nonDataChannels.size() != 1)
        {
            rdfWS_utility::messageERROR(
                "plotHists",
                "plotMode=maxVariationLines requires maxVariationChannel when more than one MC dataset is present."
            );
        }
        config.channel = nonDataChannels.front();
    }

    auto isDDIter = isDD.find(config.channel);
    if (isDDIter != isDD.end() && isDDIter->second == 1)
    {
        rdfWS_utility::messageERROR(
            "plotHists",
            "plotMode=maxVariationLines does not support DD channel " + config.channel + "."
        );
    }

    return config;
}

std::map<std::string, int> getNeedCropMap(
    rdfWS_utility::JsonObject jsonConfig,
    const std::vector<std::string> &variables)
{
    std::map<std::string, int> needCrop;
    for (const auto &varName : variables)
    {
        needCrop.emplace(varName, 0);
    }
    if (!jsonConfig.contains("needCrop"))
        return needCrop;

    auto rawConfig = jsonConfig.at("needCrop").get<nlohmann::json>();
    if (rawConfig.is_array())
    {
        for (const auto &varName : rawConfig.get<std::vector<std::string>>())
        {
            needCrop[varName] = 1;
        }
        return needCrop;
    }

    if (rawConfig.is_object())
    {
        for (auto it = rawConfig.begin(); it != rawConfig.end(); ++it)
        {
            if (!it.value().is_boolean() && !it.value().is_number_integer())
                rdfWS_utility::messageERROR(
                    "plotHists",
                    "needCrop object entry " + it.key() + " must be an integer or boolean flag."
                );
            needCrop[it.key()] = it.value().get<int>();
        }
        return needCrop;
    }

    rdfWS_utility::messageERROR(
        "plotHists",
        "needCrop must be a list of variable names or a map of variable flags, but got type "
        + std::string(rawConfig.type_name()) + " with value " + summarizeJsonValue(rawConfig)
    );

    return needCrop;
}

TH1D *cloneHistogramDetached(TH1D *source, const std::string &name)
{
    if (source == nullptr)
        return nullptr;
    auto *cloned = static_cast<TH1D *>(source->Clone(name.c_str()));
    cloned->SetDirectory(0);
    return cloned;
}

void updateEnvelopeWithCandidate(TH1D *nominal, TH1D *candidate, TH1D *upEnvelope, TH1D *downEnvelope)
{
    if (nominal == nullptr || candidate == nullptr || upEnvelope == nullptr || downEnvelope == nullptr)
        return;

    const int firstBin = 0;
    const int lastBin = nominal->GetNbinsX() + 1;
    for (int bin = firstBin; bin <= lastBin; ++bin)
    {
        const double nominalValue = nominal->GetBinContent(bin);
        const double candidateValue = candidate->GetBinContent(bin);
        const double candidateError = candidate->GetBinError(bin);

        const double currentUpDelta = upEnvelope->GetBinContent(bin) - nominalValue;
        const double currentDownDelta = downEnvelope->GetBinContent(bin) - nominalValue;
        const double candidateDelta = candidateValue - nominalValue;

        if (candidateDelta > currentUpDelta)
        {
            upEnvelope->SetBinContent(bin, candidateValue);
            upEnvelope->SetBinError(bin, candidateError);
        }
        if (candidateDelta < currentDownDelta)
        {
            downEnvelope->SetBinContent(bin, candidateValue);
            downEnvelope->SetBinError(bin, candidateError);
        }
    }
}

std::map<std::string, int> getColorScheme(
    rdfWS_utility::JsonObject jsonConfig)
{
    std::map<std::string, int> colorScheme;
    if (jsonConfig.contains("colorMapping"))
    {
        auto rawConfig = jsonConfig.at("colorMapping").get<nlohmann::json>();
        if (!rawConfig.is_object())
            rdfWS_utility::messageERROR(
                "plotHists",
                "colorMapping must be a map when provided, but got type "
                + std::string(rawConfig.type_name()) + " with value " + summarizeJsonValue(rawConfig)
            );
        colorScheme = rawConfig.get<std::map<std::string, int>>();
    }

    rdfWS_utility::JsonObject defaultConfig(
        rdfWS_utility::readJson("plotHists", "json/general_config/dataset_color.json"),
        "Default Dataset Color Config"
    );
    auto defaultColors = defaultConfig.get<std::map<std::string, int>>();
    for (const auto &[dataset, color] : defaultColors)
    {
        if (colorScheme.find(dataset) == colorScheme.end())
            colorScheme.emplace(dataset, color);
    }
    return colorScheme;
}

std::map<std::string, std::string> getChannelLabels(
    rdfWS_utility::JsonObject jsonConfig)
{
    std::map<std::string, std::string> channelLabels;
    if (jsonConfig.contains("datasetLabel"))
    {
        auto rawConfig = jsonConfig.at("datasetLabel").get<nlohmann::json>();
        if (!rawConfig.is_object())
            rdfWS_utility::messageERROR(
                "plotHists",
                "datasetLabel must be a map when provided, but got type "
                + std::string(rawConfig.type_name()) + " with value " + summarizeJsonValue(rawConfig)
            );
        channelLabels = rawConfig.get<std::map<std::string, std::string>>();
    }

    rdfWS_utility::JsonObject defaultConfig(
        rdfWS_utility::readJson("plotHists", "json/general_config/legend_label.json"),
        "Default Legend Label Config"
    );
    auto defaultLabels = defaultConfig.get<std::map<std::string, std::string>>();
    for (const auto &[dataset, label] : defaultLabels)
    {
        if (channelLabels.find(dataset) == channelLabels.end())
            channelLabels.emplace(dataset, label);
    }
    return channelLabels;
}

namespace
{
    TH1D *buildZeroHistogram(const std::string &histName, rdfWS_utility::JsonObject varConfig)
    {
        TH1D *hist = nullptr;
        if (varConfig.contains("binning"))
        {
            const auto binEdges = varConfig.at("binning").get<std::vector<double>>();
            hist = new TH1D(histName.c_str(), histName.c_str(), static_cast<int>(binEdges.size()) - 1, binEdges.data());
        }
        else
        {
            const int nBins = std::stoi(varConfig.at("nBins").get<std::string>());
            const double xMin = std::stod(varConfig.at("min").get<std::string>());
            const double xMax = std::stod(varConfig.at("max").get<std::string>());
            hist = new TH1D(histName.c_str(), histName.c_str(), nBins, xMin, xMax);
        }
        hist->SetDirectory(0);
        hist->Sumw2();
        return hist;
    }

    bool isMissingHistogramError(const std::runtime_error &ex)
    {
        const std::string message = ex.what();
        return message.find("Failed to open hist TFile:") != std::string::npos ||
               message.find(" not found in file ") != std::string::npos;
    }

    void loadHistogramOrFallback(
        HistControl &histLoader,
        const std::string &fileName,
        const std::string &histName,
        const std::string &histKey,
        float scaling,
        const std::string &varName,
        rdfWS_utility::JsonObject varConfig,
        bool allowMissingAsZero,
        bool suppressMissingWarning = false,
        const std::string &fallbackFileName = "",
        const std::string &fallbackHistName = "")
    {
        try
        {
            histLoader.loadHistogram(fileName, histName, histKey, scaling, varName);
        }
        catch (const std::runtime_error &ex)
        {
            const bool allowFallbackForMissing = allowMissingAsZero || suppressMissingWarning;
            if (!allowFallbackForMissing || !isMissingHistogramError(ex))
            {
                throw;
            }

            if (!fallbackFileName.empty() && !fallbackHistName.empty())
            {
                try
                {
                    if (!suppressMissingWarning)
                    {
                        rdfWS_utility::messageWARN(
                            "plotHists",
                            "Missing histogram input for " + histKey + " (" + fileName + ", " + histName
                            + "); using nominal fallback (" + fallbackFileName + ", " + fallbackHistName + ")."
                        );
                    }
                    histLoader.loadHistogram(fallbackFileName, fallbackHistName, histKey, scaling, varName);
                    return;
                }
                catch (const std::runtime_error &fallbackEx)
                {
                    if (!isMissingHistogramError(fallbackEx))
                        throw;
                }
            }

            if (!suppressMissingWarning)
            {
                rdfWS_utility::messageWARN("plotHists", "Missing histogram input for " + histKey + " (" + fileName + ", " + histName + "); using an all-zero placeholder.");
            }
            std::unique_ptr<TH1D> zeroHist(buildZeroHistogram(histName, varConfig));
            histLoader.addHistogram(zeroHist.get(), histKey, varName);
        }
    }

    TH1D *sumStackHists(
        HistControl &histLoader,
        const std::vector<std::string> &stackOrder,
        const std::string &histName)
    {
        if (stackOrder.empty())
            return nullptr;

        auto stackHists = histLoader.getHists(stackOrder);
        TH1D *totalHist = nullptr;
        for (const auto &name : stackOrder)
        {
            auto histIter = stackHists.find(name);
            if (histIter == stackHists.end() || histIter->second == nullptr)
                continue;
            if (totalHist == nullptr)
            {
                totalHist = static_cast<TH1D *>(histIter->second->Clone(histName.c_str()));
                totalHist->SetDirectory(0);
                totalHist->Reset();
            }
            totalHist->Add(histIter->second);
        }
        for (auto &[name, hist] : stackHists)
        {
            delete hist;
        }
        return totalHist;
    }

    TH1D *loadStackVariationTotal(
        const std::string &inDir,
        const std::vector<std::string> &campaigns,
        const std::vector<std::string> &loadChannels,
        const std::map<std::string, int> &isData,
        const std::map<std::string, double> &processScale,
        const std::map<std::string, std::vector<std::string>> &needMerge,
        const std::vector<std::string> &stackOrder,
        const std::string &varName,
        const std::string &weightName,
        const std::string &nominalWeightName,
        float mcScaling,
        rdfWS_utility::JsonObject varConfig,
        bool allowMissingAsZero,
        bool suppressMissingWarning,
        bool doCrop,
        const std::vector<double> &cropRange,
        bool doNormalize)
    {
        if (stackOrder.empty())
            return nullptr;

        std::set<std::string> channelsToLoad;
        for (const auto &stackName : stackOrder)
        {
            auto mergeIt = needMerge.find(stackName);
            if (mergeIt == needMerge.end())
            {
                channelsToLoad.insert(stackName);
            }
            else
            {
                for (const auto &component : mergeIt->second)
                    channelsToLoad.insert(component);
            }
        }

        HistControl histLoader;
        for (const auto &campaign : campaigns)
        {
            HistControl campaignLoader;
            int histNum = 0;
            for (const auto &ch : loadChannels)
            {
                if (channelsToLoad.find(ch) == channelsToLoad.end())
                    continue;
                auto isDataIter = isData.find(ch);
                if (isDataIter != isData.end() && isDataIter->second == 1)
                    continue;
                const std::string histName = ch + "_" + varName + "_" + weightName;
                const std::string fallbackHistName = ch + "_" + varName + "_" + nominalWeightName;
                loadHistogramOrFallback(
                    campaignLoader,
                    inDir + "_" + campaign + "/" + histName + ".root",
                    histName,
                    ch,
                    mcScaling * processScale.at(ch),
                    varName,
                    varConfig,
                    allowMissingAsZero,
                    suppressMissingWarning,
                    inDir + "_" + campaign + "/" + fallbackHistName + ".root",
                    fallbackHistName
                );
                histNum++;
            }
            if (histNum > 0)
                histLoader.absorbHistograms(campaignLoader);
        }

        if (doCrop)
        {
            histLoader = histLoader.cropHistograms(cropRange[0], cropRange[1]);
        }

        for (const auto &[mergedCh, components] : needMerge)
        {
            histLoader.mergeHistograms(components, mergedCh);
        }

        TH1D *totalHist = sumStackHists(histLoader, stackOrder, "stack_" + varName + "_" + weightName);
        if (totalHist != nullptr && doNormalize)
        {
            const double integral = totalHist->Integral();
            if (integral > 0.0)
                totalHist->Scale(1.0 / integral);
        }
        return totalHist;
    }

    TH1D *loadChannelVariationTotal(
        const std::string &inDir,
        const std::vector<std::string> &campaigns,
        const std::map<std::string, int> &isData,
        const std::map<std::string, double> &processScale,
        const std::map<std::string, std::vector<std::string>> &needMerge,
        const std::string &channel,
        const std::string &varName,
        const std::string &weightName,
        const std::string &nominalWeightName,
        float mcScaling,
        rdfWS_utility::JsonObject varConfig,
        bool allowMissingAsZero,
        bool suppressMissingWarning,
        bool doCrop,
        const std::vector<double> &cropRange)
    {
        std::vector<std::string> channelsToLoad = {channel};
        auto mergeIter = needMerge.find(channel);
        if (mergeIter != needMerge.end())
            channelsToLoad = mergeIter->second;

        HistControl histLoader;
        for (const auto &campaign : campaigns)
        {
            HistControl campaignLoader;
            int histNum = 0;
            for (const auto &ch : channelsToLoad)
            {
                auto isDataIter = isData.find(ch);
                if (isDataIter != isData.end() && isDataIter->second == 1)
                    continue;
                const std::string histName = ch + "_" + varName + "_" + weightName;
                const std::string fallbackHistName = ch + "_" + varName + "_" + nominalWeightName;
                loadHistogramOrFallback(
                    campaignLoader,
                    inDir + "_" + campaign + "/" + histName + ".root",
                    histName,
                    ch,
                    mcScaling * processScale.at(ch),
                    varName,
                    varConfig,
                    allowMissingAsZero,
                    suppressMissingWarning,
                    inDir + "_" + campaign + "/" + fallbackHistName + ".root",
                    fallbackHistName
                );
                histNum++;
            }
            if (histNum > 0)
                histLoader.absorbHistograms(campaignLoader);
        }

        if (doCrop)
            histLoader = histLoader.cropHistograms(cropRange[0], cropRange[1]);

        if (mergeIter != needMerge.end())
            histLoader.mergeHistograms(mergeIter->second, channel);

        auto channelHists = histLoader.getHists({channel});
        TH1D *outHist = nullptr;
        auto histIter = channelHists.find(channel);
        if (histIter != channelHists.end())
        {
            outHist = histIter->second;
            channelHists.erase(histIter);
        }
        for (auto &[name, hist] : channelHists)
            delete hist;
        return outHist;
    }
}

int main(int argc, char *argv[])
{
    // read in the config
    if (argc < 2)
    {
        rdfWS_utility::messageERROR("plotHists", "No hist plot job json provided!");
    }
    std::string jsonPath = argv[1];
    rdfWS_utility::JsonObject jsonConfig(rdfWS_utility::readJson("plotHists", jsonPath), "Job Config");

    // job type check, better to have to avoid confusion
    if (jsonConfig.contains("jobType"))
    {
        std::string jobType = jsonConfig.at("jobType").get<std::string>();
        if (jobType != "plot") 
        {
            rdfWS_utility::messageERROR("plotHists.cxx", "The jobType of your config is not plot! Please check again. Running ceases...");
            exit(1);
        }
    }

    // parse basic info from json
    std::string jobName = jsonConfig.at("name");
    std::vector<std::string> runEra = jsonConfig.at("era");
    // split MC and data era to ensure not bounded
    std::vector<std::string> mcCampaign = jsonConfig.at("mc_era");

    // get lumi-value of the plot
    std::string lumiPath = jsonConfig.at("lumiConfig");
    rdfWS_utility::JsonObject lumiConfig(rdfWS_utility::readJson("plotHists", lumiPath), "Lumi Config");
    float lumiValue = sumLumiValues(lumiConfig, runEra, "data era");

    // formulate draw texts
    // moving the CMS text and luminosity to the header
    std::vector<std::string> drawHeader = jsonConfig.at("header");
    for (int i = 0; i < drawHeader.size(); i++)
    {
        std::string text = drawHeader[i];
        if (text.find("%1.0f") != std::string::npos)
            drawHeader[i] = Form(text.c_str(), lumiValue);
    }
    std::vector<std::string> drawText = jsonConfig.at("texts");


    // MC lumi used for rescaling MC
    float mclumi = sumLumiValues(lumiConfig, mcCampaign, "mc_era");
    auto mcScaling = lumiValue / mclumi;
    std::cout << "mc lumi scaling value:" << mcScaling << std::endl;

    std::vector<std::string> signalCampaign = mcCampaign;
    float signalScaling = mcScaling;
    if (jsonConfig.contains("signal_era"))
    {
        signalCampaign = jsonConfig.at("signal_era").get<std::vector<std::string>>();
        const float signalLumi = sumLumiValues(lumiConfig, signalCampaign, "signal_era");
        signalScaling = lumiValue / signalLumi;
        std::cout << "signal lumi scaling value:" << signalScaling << std::endl;
    }

    // load channels
    std::vector<std::string> channels = jsonConfig.at("datasets");
    std::map<std::string, std::string> channelLabels = getChannelLabels(jsonConfig);
    std::map<std::string, std::vector<std::string>> needMerge = jsonConfig.at("needMerge");
    std::vector<std::string> loadChannels;
    for (const auto &ch : channels)
    {
        if (needMerge.find(ch) == needMerge.end())
        {
            loadChannels.push_back(ch);
        }
        else
        {
            auto mergeList = needMerge[ch];
            for (auto mergeCh : mergeList)
            {
                loadChannels.push_back(mergeCh);
            }
        }
    }

    // for default stack plots
    std::vector<std::string> stackOrder = jsonConfig.at("stackOrder");
    int reOrder = jsonConfig.at("reOrder");
    std::vector<std::string> numerator = jsonConfig.at("numerator");

    // treat styles
    const std::vector<std::string> isSignal = jsonConfig.at("isSignal");
    const std::set<std::string> signalLoadChannels = expandSignalLoadChannels(isSignal, needMerge);
    std::map<std::string, int> isData = getIsDataMap(jsonConfig, channels);
    const std::map<std::string, int> isDD = getIsDDMap(jsonConfig, channels);
    const std::map<std::string, double> ddVariations = getDDVariationMap(jsonConfig, channels);
    const std::map<std::string, double> processScale = getProcessScaleMap(jsonConfig, loadChannels);
    bool withDD = false;
    for (const auto &[channel, flag] : isDD)
    {
        if (flag == 1)
        {
            withDD = true;
            break;
        }
    }
    const OverlaySystConfig overlaySystConfig = getOverlaySystConfig(jsonConfig, channels, stackOrder, isData, isDD);
    std::string dataWeight = jsonConfig.at("dataWeight");
    std::string MCWeight = jsonConfig.at("MCWeight");
    std::string DDWeight = "";
    if (jsonConfig.contains("DDWeight"))
    {
        std::string tempWeight = jsonConfig.at("DDWeight");
        DDWeight = tempWeight;
    }
    std::vector<std::string> ddCampaign;
    if (withDD)
    {
        if (!jsonConfig.contains("dd_era"))
            rdfWS_utility::messageERROR("plotHists", "isDD is enabled but dd_era is not provided.");
        if (DDWeight == "")
            rdfWS_utility::messageERROR("plotHists", "isDD is enabled but DDWeight is not provided.");
        ddCampaign = jsonConfig.at("dd_era").get<std::vector<std::string>>();
    }
    const MaxVariationLinePlotConfig maxVariationLinePlot = getMaxVariationLinePlotConfig(jsonConfig, channels, isData, isDD);
    float ddScaling(0.0);
    if (withDD)
    {
        float ddlumi = sumLumiValues(lumiConfig, ddCampaign, "dd_era");
        ddScaling = lumiValue / ddlumi;
        std::cout << "data-driven lumi scaling value:" << ddScaling << std::endl;
    }
    const std::map<std::string, int> colorScheme = getColorScheme(jsonConfig);

    std::string inDir = jsonConfig.at("inDir");
    bool allowMissingAsZero = false;
    if (jsonConfig.contains("allowMissingAsZero"))
    {
        allowMissingAsZero = jsonConfig.at("allowMissingAsZero");
    }
    bool silentMissingSystForCrossEra = false;
    if (jsonConfig.contains("silentMissingSystForCrossEra"))
    {
        silentMissingSystForCrossEra = jsonConfig.at("silentMissingSystForCrossEra");
    }
    bool drawStackSystBand = true;
    if (jsonConfig.contains("drawStackSystBand"))
    {
        drawStackSystBand = jsonConfig.at("drawStackSystBand").get<bool>();
    }

    // for each variable, load the histograms and plot
    std::vector<std::string> variables = jsonConfig.at("varNames");
    std::map<std::string, int> needCrop = getNeedCropMap(jsonConfig, variables);
    for (const auto &varName : variables)
    {
        // setup drawing options
        PlotContext options;
        options.doLog = jsonConfig.at("doLogPlot");
        options.doNormalize = jsonConfig.at("doNormPlot");

        std::vector<std::string> isSignalCopy = isSignal;
        options.isSignal = isSignalCopy;
        std::map<std::string, int> isDataCopy = isData;
        std::map<std::string, int> isDDCopy = isDD;
        options.isData.push_back(isDataCopy);
        std::string varConfigPath = jsonConfig.at("varConfig");
        rdfWS_utility::JsonObject varJson(rdfWS_utility::readJson("plotHists", varConfigPath), "Var Config");
        auto varConfig = varJson.at(varName);
        options.xLabel = varConfig.at("label").get<std::string>();
        options.displayEqualWidthBins = jsonConfig.contains("displayEqualWidthBins")
                                           ? jsonConfig.at("displayEqualWidthBins").get<bool>()
                                           : false;
        if (options.displayEqualWidthBins)
        {
            if (!varConfig.contains("binning"))
            {
                rdfWS_utility::messageERROR("plotHists", "displayEqualWidthBins requires varConfig.binning for variable " + varName + ".");
            }
            options.displayBinEdges = varConfig.at("binning").get<std::vector<double>>();
        }
        options.blindRanges = getBlindRanges(jsonConfig, varName);

        // add default values for y-axis labels
        std::string yUpperLabel = "Events";
        if (options.doNormalize) yUpperLabel = "a.u.";
        std::string yLowerLabel = "Data / MC";
        if (jsonConfig.contains("yLabel"))
        {
            yUpperLabel = jsonConfig.at("yLabel").get<std::string>();
        }
        if (jsonConfig.contains("yRatioLabe"))
        {
            yLowerLabel = jsonConfig.at("yRatioLabel").get<std::string>();
        }
        options.yLabel = {yUpperLabel, yLowerLabel};
        options.xSize = jsonConfig.at("histXSize");
        options.ySize = jsonConfig.at("histYSize");

        std::vector<std::string> binLabels = varConfig.at("binLabels");
        if (options.displayEqualWidthBins)
        {
            binLabels.clear();
        }
        // need to truncate if have cropping
        if (needCrop[varName] == 1 && binLabels.size()>0)
        {
            std::vector<double> cropRange = jsonConfig.at("cropedRange").at(varName);
            const auto selectedBins = getCropBinIndicesFromConfig(varConfig, cropRange);
            if (selectedBins.empty())
                rdfWS_utility::messageERROR("plotHists", "Crop range does not fully contain any bin for variable " + varName + ".");

            std::vector<std::string> croppedLabels;
            croppedLabels.reserve(selectedBins.size());
            for (const int index : selectedBins)
            {
                if (index < 0 || index >= static_cast<int>(binLabels.size()))
                    rdfWS_utility::messageERROR("plotHists", "Cropped bin-label index out of range for variable " + varName + ".");
                croppedLabels.push_back(binLabels[index]);
            }
            binLabels = croppedLabels;
        }

        // configure output
        std::string outputDir = jsonConfig.at("outDir");
        rdfWS_utility::creatingFolder("plotHists", outputDir);
        std::string plotName = outputDir + "/";
        if (options.doNormalize)
            plotName += "norm_";
        plotName += "data_MC";
        for (auto era : runEra)
        {
            plotName += "_";
            plotName += era;
        }
        plotName += "_" + varName + "_" + jobName;
        if (needCrop[varName] == 1)
            plotName += "_crop";
        PlotControl pHelper(plotName);

        // read the histograms separatedly with MC and data
        HistControl histLoader;
        // first data era: runEra[0]
        for (const auto &ch : loadChannels)
        {
            std::string histName = ch + "_" + varName + "_";
            if (isDataCopy[ch] == 1)
                histName += dataWeight;
            else
                continue;
            loadHistogramOrFallback(histLoader, inDir + "_" + runEra[0] + "/" + histName + ".root", histName, ch, 1.0, varName, varConfig, allowMissingAsZero);
        }

        HistControl mcHistLoader;
        int withMC=0;
        for (const auto &ch : loadChannels)
        {
            if (isDataCopy[ch] == 1)
                continue;
            if (isDDCopy[ch] == 1)
                continue;
            const bool isSignalLoadChannel = signalLoadChannels.find(ch) != signalLoadChannels.end();
            const auto &campaignsToLoad = isSignalLoadChannel ? signalCampaign : mcCampaign;
            const float scalingToUse = isSignalLoadChannel ? signalScaling : mcScaling;
            HistControl channelLoader;
            for (const auto &campaign : campaignsToLoad)
            {
                HistControl campaignLoader;
                std::string histName = ch + "_" + varName + "_" + MCWeight;
                loadHistogramOrFallback(campaignLoader, inDir + "_" + campaign + "/" + histName + ".root", histName, ch, scalingToUse * processScale.at(ch), varName, varConfig, allowMissingAsZero);
                channelLoader.absorbHistograms(campaignLoader);
            }
            mcHistLoader.absorbHistograms(channelLoader);
            withMC++;
        }
        if (withMC>0)
            histLoader.absorbHistograms(mcHistLoader);

        if (withDD)
        {
            HistControl DDHistLoader;
            for (const auto &ch : loadChannels)
            {
                if (isDDCopy[ch] != 1)
                    continue;
                std::string histName = ch + "_" + varName + "_" + DDWeight;
                loadHistogramOrFallback(DDHistLoader, inDir + "_" + ddCampaign[0] + "/" + histName + ".root", histName, ch, ddScaling * processScale.at(ch), varName, varConfig, allowMissingAsZero);
            }
            histLoader.absorbHistograms(DDHistLoader);
        }

        // add together when more than 1 eras
        for (int i = 1; i < runEra.size(); i++)
        {
            HistControl tempHistLoader;
            int histNum=0;
            for (const auto &ch : loadChannels)
            {
                std::string histName = ch + "_" + varName + "_";
                if (isDataCopy[ch] == 1)
                    histName += dataWeight;
                else if (isDDCopy[ch] == 1)
                    continue;
                else
                    continue;
                loadHistogramOrFallback(tempHistLoader, inDir + "_" + runEra[i] + "/" + histName + ".root", histName, ch, 1.0, varName, varConfig, allowMissingAsZero);
                histNum++;
            }
            if (histNum>0) 
                histLoader.absorbHistograms(tempHistLoader);
        }

        // Nov 26, add DD contribution
        if (withDD)
        {
            for (int i = 1; i < ddCampaign.size(); i++)
            {
                HistControl tempDDHistLoader;
                int histNum=0;
                for (const auto &ch : loadChannels)
                {
                    if (isDDCopy[ch] != 1)
                        continue;
                    std::string histName = ch + "_" + varName + "_" + DDWeight;
                    loadHistogramOrFallback(tempDDHistLoader, inDir + "_" + ddCampaign[i] + "/" + histName + ".root", histName, ch, ddScaling * processScale.at(ch), varName, varConfig, allowMissingAsZero);
                    histNum++;
                }
                if (histNum>0)
                    histLoader.absorbHistograms(tempDDHistLoader);
            }
        }

        // take crops TODO
        if (needCrop[varName] == 1)
        {
            auto loadedHists = histLoader.getHistInstance();
            if (loadedHists.empty())
            {
                rdfWS_utility::messageERROR("plotHists", "No histograms loaded for crop checking.");
            }
            auto refHist = loadedHists.begin()->second;
            std::vector<double> cropRange = jsonConfig.at("cropedRange").at(varName);
            const auto selectedBins = getCropBinIndicesFromConfig(varConfig, cropRange);
            if (selectedBins.empty())
                rdfWS_utility::messageERROR("plotHists", "Crop range does not fully contain any bin for variable " + varName + ".");

            if (varConfig.contains("binning"))
            {
                const auto cfgEdges = varConfig.at("binning").get<std::vector<double>>();
                if (refHist->GetNbinsX() != static_cast<int>(cfgEdges.size()) - 1)
                    rdfWS_utility::messageERROR("plotHists", "Histogram binning and varConfig disagree for cropped variable " + varName + ".");
                for (int i = 1; i <= refHist->GetNbinsX(); ++i)
                {
                    if (std::abs(refHist->GetXaxis()->GetBinLowEdge(i) - cfgEdges[i - 1]) > 1e-9 ||
                        std::abs(refHist->GetXaxis()->GetBinUpEdge(i) - cfgEdges[i]) > 1e-9)
                    {
                        rdfWS_utility::messageERROR("plotHists", "Histogram binning and varConfig disagree for cropped variable " + varName + ".");
                    }
                }
            }
            else
            {
                int cfgNBins = std::stoi(varConfig.at("nBins").get<std::string>());
                double cfgMin = std::stod(varConfig.at("min").get<std::string>());
                double cfgMax = std::stod(varConfig.at("max").get<std::string>());
                if (refHist->GetNbinsX() != cfgNBins || std::abs(refHist->GetXaxis()->GetXmin() - cfgMin) > 1e-9 || std::abs(refHist->GetXaxis()->GetXmax() - cfgMax) > 1e-9)
                {
                    rdfWS_utility::messageERROR("plotHists", "Histogram binning and varConfig disagree for cropped variable " + varName + ". Crop would not produce the expected result.");
                }
            }
            histLoader = histLoader.cropHistograms(cropRange[0], cropRange[1]);
        }

        // merge the plots according to the config
        for (const auto &[mergedCh, components] : needMerge)
        {
            histLoader.mergeHistograms(components, mergedCh);
        }

        auto histsNeeded = histLoader.getHists(channels);
        int doRatio = jsonConfig.at("doRatioPlot");
        options.isData.push_back(std::map<std::string, int>{});
        for (auto key : numerator)
        {
            options.isData[1].emplace(key, isDataCopy[key]);
        }

        std::map<std::string, TH1D *> ratioHists = {};
        // currently only allowing ratio plot with respect to stack
        if (doRatio && stackOrder.size() > 0)
        {
            ratioHists = histLoader.getRatios(numerator, stackOrder, options.doNormalize);
        }

        // to add possible uncertainties
        auto stackUncert = drawStackSystBand ? getStackUncert(stackOrder, jsonConfig) : std::vector<std::map<std::string, double>>{};
        std::map<std::string, double> stackUp, stackDown;
        if (stackUncert.size() == 2)
        {
            stackUp = stackUncert[0];
            stackDown = stackUncert[1];
        }
        std::vector<std::pair<TH1D *, TH1D *>> shapeUncerts;
        const auto weightSysts = getWeightSystematics(jsonConfig);
        const auto signalWeightSysts = jsonConfig.contains("signalWeightSyst")
            ? getWeightSystematics(jsonConfig, "signalWeightSyst")
            : weightSysts;
        const auto ddWeightSysts = jsonConfig.contains("DDweightSyst")
            ? getWeightSystematics(jsonConfig, "DDweightSyst")
            : std::map<std::string, std::vector<std::string>>{};
        const auto shiftSystNames = getShiftSystematicNames(jsonConfig);
        const auto signalShiftSystNames =
            (jsonConfig.contains("signalShiftSyst") || jsonConfig.contains("signalShiftExprSyst"))
                ? getShiftSystematicNames(jsonConfig, "signalShiftSyst", "signalShiftExprSyst")
                : shiftSystNames;
        const auto ddShiftSystNames =
            (jsonConfig.contains("DDshiftSyst") || jsonConfig.contains("DDshiftExprSyst"))
                ? getShiftSystematicNames(jsonConfig, "DDshiftSyst", "DDshiftExprSyst")
                : std::vector<std::string>{};

        if (maxVariationLinePlot.enabled)
        {
            auto nominalIter = histsNeeded.find(maxVariationLinePlot.channel);
            if (nominalIter == histsNeeded.end() || nominalIter->second == nullptr)
            {
                rdfWS_utility::messageERROR(
                    "plotHists",
                    "plotMode=maxVariationLines could not find nominal histogram for channel " + maxVariationLinePlot.channel + "."
                );
            }

            const bool isSignalChannel = std::find(isSignal.begin(), isSignal.end(), maxVariationLinePlot.channel) != isSignal.end();
            const bool isDDChannel = (isDDCopy[maxVariationLinePlot.channel] == 1);
            const auto &campaignsToLoad = isSignalChannel ? signalCampaign : (isDDChannel ? ddCampaign : mcCampaign);
            const float scalingToUse = isSignalChannel ? signalScaling : (isDDChannel ? ddScaling : mcScaling);
            const std::string nominalWeightToUse = isSignalChannel ? MCWeight : (isDDChannel ? DDWeight : MCWeight);
            const auto &channelWeightSysts = isSignalChannel ? signalWeightSysts : (isDDChannel ? ddWeightSysts : weightSysts);
            const auto &channelShiftSystNames = isSignalChannel ? signalShiftSystNames : (isDDChannel ? ddShiftSystNames : shiftSystNames);

            std::vector<double> cropRange = {};
            if (needCrop[varName] == 1)
                cropRange = jsonConfig.at("cropedRange").at(varName).get<std::vector<double>>();

            const std::string nominalKey = maxVariationLinePlot.channel + "__nominal";
            const std::string upKey = maxVariationLinePlot.channel + "__syst_a_up";
            const std::string downKey = maxVariationLinePlot.channel + "__syst_b_down";

            TH1D *nominalHist = cloneHistogramDetached(nominalIter->second, nominalKey);
            TH1D *upEnvelope = cloneHistogramDetached(nominalIter->second, upKey);
            TH1D *downEnvelope = cloneHistogramDetached(nominalIter->second, downKey);

            for (const auto &[systName, variations] : channelWeightSysts)
            {
                const std::string upWeightName = nominalWeightToUse + "_" + systName + "_up";
                const std::string downWeightName = nominalWeightToUse + "_" + systName + "_down";
                TH1D *upHist = loadChannelVariationTotal(inDir, campaignsToLoad, isDataCopy, processScale, needMerge, maxVariationLinePlot.channel, varName, upWeightName, nominalWeightToUse, scalingToUse, varConfig, allowMissingAsZero, silentMissingSystForCrossEra, needCrop[varName] == 1, cropRange);
                TH1D *downHist = loadChannelVariationTotal(inDir, campaignsToLoad, isDataCopy, processScale, needMerge, maxVariationLinePlot.channel, varName, downWeightName, nominalWeightToUse, scalingToUse, varConfig, allowMissingAsZero, silentMissingSystForCrossEra, needCrop[varName] == 1, cropRange);
                updateEnvelopeWithCandidate(nominalHist, upHist, upEnvelope, downEnvelope);
                updateEnvelopeWithCandidate(nominalHist, downHist, upEnvelope, downEnvelope);
                delete upHist;
                delete downHist;
            }
            for (const auto &systName : channelShiftSystNames)
            {
                const std::string upWeightName = nominalWeightToUse + "_" + systName + "_up";
                const std::string downWeightName = nominalWeightToUse + "_" + systName + "_down";
                TH1D *upHist = loadChannelVariationTotal(inDir, campaignsToLoad, isDataCopy, processScale, needMerge, maxVariationLinePlot.channel, varName, upWeightName, nominalWeightToUse, scalingToUse, varConfig, allowMissingAsZero, silentMissingSystForCrossEra, needCrop[varName] == 1, cropRange);
                TH1D *downHist = loadChannelVariationTotal(inDir, campaignsToLoad, isDataCopy, processScale, needMerge, maxVariationLinePlot.channel, varName, downWeightName, nominalWeightToUse, scalingToUse, varConfig, allowMissingAsZero, silentMissingSystForCrossEra, needCrop[varName] == 1, cropRange);
                updateEnvelopeWithCandidate(nominalHist, upHist, upEnvelope, downEnvelope);
                updateEnvelopeWithCandidate(nominalHist, downHist, upEnvelope, downEnvelope);
                delete upHist;
                delete downHist;
            }

            std::map<std::string, TH1D *> lineHists;
            lineHists.emplace(nominalKey, nominalHist);
            lineHists.emplace(upKey, upEnvelope);
            lineHists.emplace(downKey, downEnvelope);

            PlotContext lineOptions = options;
            lineOptions.isData.clear();
            lineOptions.isData.push_back({
                {nominalKey, 0},
                {upKey, 0},
                {downKey, 0},
            });
            lineOptions.isSignal.clear();

            std::map<std::string, int> lineColorScheme = colorScheme;
            lineColorScheme[nominalKey] = 1;
            lineColorScheme[upKey] = 632;
            lineColorScheme[downKey] = 600;

            std::map<std::string, std::string> lineLabels = channelLabels;
            const auto labelIter = channelLabels.find(maxVariationLinePlot.channel);
            lineLabels[nominalKey] = (labelIter != channelLabels.end()) ? labelIter->second : maxVariationLinePlot.channel;
            lineLabels[upKey] = "syst up";
            lineLabels[downKey] = "syst down";

            pHelper.drawStackHistWithRatio(
                lineHists,
                {},
                {},
                {},
                0,
                {},
                lineOptions,
                1.0,
                lineColorScheme,
                lineLabels,
                drawHeader,
                drawText,
                {},
                {},
                {},
                false,
                binLabels
            );

            for (auto &[histName, hist] : lineHists)
                delete hist;
            lineHists.clear();

            for (auto &[histName, hist] : histsNeeded)
                delete hist;
            histsNeeded.clear();
            continue;
        }

        std::vector<std::string> mcWeightStackOrder;
        for (const auto &stackName : stackOrder)
        {
            bool isSignalChannel = false;
            for (const auto &signalChannel : isSignal)
            {
                if (stackName == signalChannel)
                {
                    isSignalChannel = true;
                    break;
                }
            }
            if (isDDCopy[stackName] != 1 && !isSignalChannel)
                mcWeightStackOrder.push_back(stackName);
        }

        // MC variations are loaded from mcWeightStackOrder, which deliberately
        // excludes DD channels. Add the nominal DD contribution back before
        // comparing the variation with the full nominal stack.
        auto addNominalDDToVariation = [&](TH1D *variation)
        {
            if (variation == nullptr)
                return;
            for (const auto &ddChannel : stackOrder)
            {
                if (isDDCopy[ddChannel] != 1)
                    continue;
                auto nominalDDIter = histsNeeded.find(ddChannel);
                if (nominalDDIter != histsNeeded.end() && nominalDDIter->second != nullptr)
                    variation->Add(nominalDDIter->second);
            }
        };

        if (drawStackSystBand) for (const auto &[systName, variations] : weightSysts)
        {
            if (variations.size() != 2 && variations.size() != 3)
            {
                rdfWS_utility::messageERROR("plotHists", "weightSyst entry " + systName + " must be [up, down] or [nominal, up, down].");
                exit(1);
            }
            const std::string upWeightName = MCWeight + "_" + systName + "_up";
            const std::string downWeightName = MCWeight + "_" + systName + "_down";
            std::vector<double> cropRange = {};
            if (needCrop[varName] == 1)
                cropRange = jsonConfig.at("cropedRange").at(varName).get<std::vector<double>>();
            TH1D *upHist = loadStackVariationTotal(inDir, mcCampaign, loadChannels, isDataCopy, processScale, needMerge, mcWeightStackOrder, varName, upWeightName, MCWeight, mcScaling, varConfig, allowMissingAsZero, silentMissingSystForCrossEra, needCrop[varName] == 1, cropRange, false);
            TH1D *downHist = loadStackVariationTotal(inDir, mcCampaign, loadChannels, isDataCopy, processScale, needMerge, mcWeightStackOrder, varName, downWeightName, MCWeight, mcScaling, varConfig, allowMissingAsZero, silentMissingSystForCrossEra, needCrop[varName] == 1, cropRange, false);
            addNominalDDToVariation(upHist);
            addNominalDDToVariation(downHist);
            if (options.doNormalize)
            {
                if (upHist != nullptr && upHist->Integral() > 0.0)
                    upHist->Scale(1.0 / upHist->Integral());
                if (downHist != nullptr && downHist->Integral() > 0.0)
                    downHist->Scale(1.0 / downHist->Integral());
            }
            if (upHist != nullptr && downHist != nullptr)
            {
                shapeUncerts.emplace_back(upHist, downHist);
            }
            else
            {
                delete upHist;
                delete downHist;
            }
        }
        if (drawStackSystBand) for (const auto &[systName, variations] : ddWeightSysts)
        {
            if (variations.size() != 2 && variations.size() != 3)
            {
                rdfWS_utility::messageERROR("plotHists", "DDweightSyst entry " + systName + " must be [up, down] or [nominal, up, down].");
                exit(1);
            }
            const std::string upWeightName = DDWeight + "_" + systName + "_up";
            const std::string downWeightName = DDWeight + "_" + systName + "_down";
            std::vector<double> cropRange = {};
            if (needCrop[varName] == 1)
                cropRange = jsonConfig.at("cropedRange").at(varName).get<std::vector<double>>();
            TH1D *upHist = loadStackVariationTotal(inDir, mcCampaign, loadChannels, isDataCopy, processScale, needMerge, mcWeightStackOrder, varName, MCWeight, MCWeight, mcScaling, varConfig, allowMissingAsZero, silentMissingSystForCrossEra, needCrop[varName] == 1, cropRange, false);
            TH1D *downHist = loadStackVariationTotal(inDir, mcCampaign, loadChannels, isDataCopy, processScale, needMerge, mcWeightStackOrder, varName, MCWeight, MCWeight, mcScaling, varConfig, allowMissingAsZero, silentMissingSystForCrossEra, needCrop[varName] == 1, cropRange, false);
            for (const auto &ddChannel : stackOrder)
            {
                if (isDDCopy[ddChannel] != 1)
                    continue;
                TH1D *ddUpHist = loadChannelVariationTotal(inDir, ddCampaign, isDataCopy, processScale, needMerge, ddChannel, varName, upWeightName, DDWeight, ddScaling, varConfig, allowMissingAsZero, silentMissingSystForCrossEra, needCrop[varName] == 1, cropRange);
                TH1D *ddDownHist = loadChannelVariationTotal(inDir, ddCampaign, isDataCopy, processScale, needMerge, ddChannel, varName, downWeightName, DDWeight, ddScaling, varConfig, allowMissingAsZero, silentMissingSystForCrossEra, needCrop[varName] == 1, cropRange);
                auto nominalDDIter = histsNeeded.find(ddChannel);
                if (ddUpHist == nullptr || ddDownHist == nullptr)
                {
                    delete ddUpHist;
                    delete ddDownHist;
                    if (nominalDDIter == histsNeeded.end())
                        continue;
                    if (upHist == nullptr)
                    {
                        upHist = static_cast<TH1D *>(nominalDDIter->second->Clone(("stack_" + varName + "_" + upWeightName).c_str()));
                        upHist->SetDirectory(0);
                        upHist->Reset();
                    }
                    if (downHist == nullptr)
                    {
                        downHist = static_cast<TH1D *>(nominalDDIter->second->Clone(("stack_" + varName + "_" + downWeightName).c_str()));
                        downHist->SetDirectory(0);
                        downHist->Reset();
                    }
                    upHist->Add(nominalDDIter->second);
                    downHist->Add(nominalDDIter->second);
                    continue;
                }
                upHist->Add(ddUpHist);
                downHist->Add(ddDownHist);
                delete ddUpHist;
                delete ddDownHist;
            }
            if (options.doNormalize)
            {
                if (upHist != nullptr && upHist->Integral() > 0.0)
                    upHist->Scale(1.0 / upHist->Integral());
                if (downHist != nullptr && downHist->Integral() > 0.0)
                    downHist->Scale(1.0 / downHist->Integral());
            }
            if (upHist != nullptr && downHist != nullptr)
            {
                shapeUncerts.emplace_back(upHist, downHist);
            }
            else
            {
                delete upHist;
                delete downHist;
            }
        }
        if (drawStackSystBand) for (const auto &systName : shiftSystNames)
        {
            const std::string upWeightName = MCWeight + "_" + systName + "_up";
            const std::string downWeightName = MCWeight + "_" + systName + "_down";
            std::vector<double> cropRange = {};
            if (needCrop[varName] == 1)
                cropRange = jsonConfig.at("cropedRange").at(varName).get<std::vector<double>>();
            TH1D *upHist = loadStackVariationTotal(inDir, mcCampaign, loadChannels, isDataCopy, processScale, needMerge, mcWeightStackOrder, varName, upWeightName, MCWeight, mcScaling, varConfig, allowMissingAsZero, silentMissingSystForCrossEra, needCrop[varName] == 1, cropRange, false);
            TH1D *downHist = loadStackVariationTotal(inDir, mcCampaign, loadChannels, isDataCopy, processScale, needMerge, mcWeightStackOrder, varName, downWeightName, MCWeight, mcScaling, varConfig, allowMissingAsZero, silentMissingSystForCrossEra, needCrop[varName] == 1, cropRange, false);
            addNominalDDToVariation(upHist);
            addNominalDDToVariation(downHist);
            if (options.doNormalize)
            {
                if (upHist != nullptr && upHist->Integral() > 0.0)
                    upHist->Scale(1.0 / upHist->Integral());
                if (downHist != nullptr && downHist->Integral() > 0.0)
                    downHist->Scale(1.0 / downHist->Integral());
            }
            if (upHist != nullptr && downHist != nullptr)
            {
                shapeUncerts.emplace_back(upHist, downHist);
            }
            else
            {
                delete upHist;
                delete downHist;
            }
        }
        if (drawStackSystBand) for (const auto &systName : ddShiftSystNames)
        {
            const std::string upWeightName = DDWeight + "_" + systName + "_up";
            const std::string downWeightName = DDWeight + "_" + systName + "_down";
            std::vector<double> cropRange = {};
            if (needCrop[varName] == 1)
                cropRange = jsonConfig.at("cropedRange").at(varName).get<std::vector<double>>();
            TH1D *upHist = loadStackVariationTotal(inDir, mcCampaign, loadChannels, isDataCopy, processScale, needMerge, mcWeightStackOrder, varName, MCWeight, MCWeight, mcScaling, varConfig, allowMissingAsZero, silentMissingSystForCrossEra, needCrop[varName] == 1, cropRange, false);
            TH1D *downHist = loadStackVariationTotal(inDir, mcCampaign, loadChannels, isDataCopy, processScale, needMerge, mcWeightStackOrder, varName, MCWeight, MCWeight, mcScaling, varConfig, allowMissingAsZero, silentMissingSystForCrossEra, needCrop[varName] == 1, cropRange, false);
            for (const auto &ddChannel : stackOrder)
            {
                if (isDDCopy[ddChannel] != 1)
                    continue;
                TH1D *ddUpHist = loadChannelVariationTotal(inDir, ddCampaign, isDataCopy, processScale, needMerge, ddChannel, varName, upWeightName, DDWeight, ddScaling, varConfig, allowMissingAsZero, silentMissingSystForCrossEra, needCrop[varName] == 1, cropRange);
                TH1D *ddDownHist = loadChannelVariationTotal(inDir, ddCampaign, isDataCopy, processScale, needMerge, ddChannel, varName, downWeightName, DDWeight, ddScaling, varConfig, allowMissingAsZero, silentMissingSystForCrossEra, needCrop[varName] == 1, cropRange);
                auto nominalDDIter = histsNeeded.find(ddChannel);
                if (ddUpHist == nullptr || ddDownHist == nullptr)
                {
                    delete ddUpHist;
                    delete ddDownHist;
                    if (nominalDDIter == histsNeeded.end())
                        continue;
                    if (upHist == nullptr)
                    {
                        upHist = static_cast<TH1D *>(nominalDDIter->second->Clone(("stack_" + varName + "_" + upWeightName).c_str()));
                        upHist->SetDirectory(0);
                        upHist->Reset();
                    }
                    if (downHist == nullptr)
                    {
                        downHist = static_cast<TH1D *>(nominalDDIter->second->Clone(("stack_" + varName + "_" + downWeightName).c_str()));
                        downHist->SetDirectory(0);
                        downHist->Reset();
                    }
                    upHist->Add(nominalDDIter->second);
                    downHist->Add(nominalDDIter->second);
                    continue;
                }
                upHist->Add(ddUpHist);
                downHist->Add(ddDownHist);
                delete ddUpHist;
                delete ddDownHist;
            }
            if (options.doNormalize)
            {
                if (upHist != nullptr && upHist->Integral() > 0.0)
                    upHist->Scale(1.0 / upHist->Integral());
                if (downHist != nullptr && downHist->Integral() > 0.0)
                    downHist->Scale(1.0 / downHist->Integral());
            }
            if (upHist != nullptr && downHist != nullptr)
            {
                shapeUncerts.emplace_back(upHist, downHist);
            }
            else
            {
                delete upHist;
                delete downHist;
            }
        }
        if (drawStackSystBand)
        {
            for (const auto &ddChannel : stackOrder)
            {
                if (isDDCopy[ddChannel] != 1)
                    continue;

                auto ddVarIter = ddVariations.find(ddChannel);
                if (ddVarIter == ddVariations.end())
                    continue;
                const double relVariation = ddVarIter->second;
                if (relVariation <= 0.0)
                    continue;

                TH1D *upHist = nullptr;
                TH1D *downHist = nullptr;
                for (const auto &stackName : stackOrder)
                {
                    auto nominalIter = histsNeeded.find(stackName);
                    if (nominalIter == histsNeeded.end())
                        continue;

                    if (upHist == nullptr)
                    {
                        upHist = static_cast<TH1D *>(
                            nominalIter->second->Clone(("stack_" + varName + "_" + ddChannel + "_DDvariation_up").c_str())
                        );
                        upHist->SetDirectory(0);
                        upHist->Reset();
                    }
                    if (downHist == nullptr)
                    {
                        downHist = static_cast<TH1D *>(
                            nominalIter->second->Clone(("stack_" + varName + "_" + ddChannel + "_DDvariation_down").c_str())
                        );
                        downHist->SetDirectory(0);
                        downHist->Reset();
                    }

                    const double upScale = (stackName == ddChannel) ? (1.0 + relVariation) : 1.0;
                    const double downScale = (stackName == ddChannel) ? std::max(0.0, 1.0 - relVariation) : 1.0;
                    upHist->Add(nominalIter->second, upScale);
                    downHist->Add(nominalIter->second, downScale);
                }

                if (options.doNormalize)
                {
                    if (upHist != nullptr && upHist->Integral() > 0.0)
                        upHist->Scale(1.0 / upHist->Integral());
                    if (downHist != nullptr && downHist->Integral() > 0.0)
                        downHist->Scale(1.0 / downHist->Integral());
                }

                if (upHist != nullptr && downHist != nullptr)
                    shapeUncerts.emplace_back(upHist, downHist);
                else
                {
                    delete upHist;
                    delete downHist;
                }
            }
        }
        std::map<std::string, std::vector<std::pair<TH1D *, TH1D *>>> nonStackShapeUncerts;
        for (const auto &channel : overlaySystConfig.channels)
        {
            auto isDataIter = isDataCopy.find(channel);
            if (isDataIter != isDataCopy.end() && isDataIter->second == 1)
                continue;
            const bool isSignalChannel = std::find(isSignal.begin(), isSignal.end(), channel) != isSignal.end();
            const bool isDDChannel = (isDDCopy[channel] == 1);
            const auto &campaignsToLoad = isSignalChannel ? signalCampaign : (isDDChannel ? ddCampaign : mcCampaign);
            const float scalingToUse = isSignalChannel ? signalScaling : (isDDChannel ? ddScaling : mcScaling);
            const std::string nominalWeightToUse = isSignalChannel ? MCWeight : (isDDChannel ? DDWeight : MCWeight);

            std::vector<std::pair<TH1D *, TH1D *>> channelShapeUncerts;
            std::vector<double> cropRange = {};
            if (needCrop[varName] == 1)
                cropRange = jsonConfig.at("cropedRange").at(varName).get<std::vector<double>>();

            const auto &channelWeightSysts = isSignalChannel ? signalWeightSysts : (isDDChannel ? ddWeightSysts : weightSysts);
            const auto &channelShiftSystNames = isSignalChannel ? signalShiftSystNames : (isDDChannel ? ddShiftSystNames : shiftSystNames);

            for (const auto &[systName, variations] : channelWeightSysts)
            {
                const std::string upWeightName = nominalWeightToUse + "_" + systName + "_up";
                const std::string downWeightName = nominalWeightToUse + "_" + systName + "_down";
                TH1D *upHist = loadChannelVariationTotal(inDir, campaignsToLoad, isDataCopy, processScale, needMerge, channel, varName, upWeightName, nominalWeightToUse, scalingToUse, varConfig, allowMissingAsZero, silentMissingSystForCrossEra, needCrop[varName] == 1, cropRange);
                TH1D *downHist = loadChannelVariationTotal(inDir, campaignsToLoad, isDataCopy, processScale, needMerge, channel, varName, downWeightName, nominalWeightToUse, scalingToUse, varConfig, allowMissingAsZero, silentMissingSystForCrossEra, needCrop[varName] == 1, cropRange);
                if (upHist != nullptr && downHist != nullptr)
                    channelShapeUncerts.emplace_back(upHist, downHist);
                else
                {
                    delete upHist;
                    delete downHist;
                }
            }
            for (const auto &systName : channelShiftSystNames)
            {
                const std::string upWeightName = nominalWeightToUse + "_" + systName + "_up";
                const std::string downWeightName = nominalWeightToUse + "_" + systName + "_down";
                TH1D *upHist = loadChannelVariationTotal(inDir, campaignsToLoad, isDataCopy, processScale, needMerge, channel, varName, upWeightName, nominalWeightToUse, scalingToUse, varConfig, allowMissingAsZero, silentMissingSystForCrossEra, needCrop[varName] == 1, cropRange);
                TH1D *downHist = loadChannelVariationTotal(inDir, campaignsToLoad, isDataCopy, processScale, needMerge, channel, varName, downWeightName, nominalWeightToUse, scalingToUse, varConfig, allowMissingAsZero, silentMissingSystForCrossEra, needCrop[varName] == 1, cropRange);
                if (upHist != nullptr && downHist != nullptr)
                    channelShapeUncerts.emplace_back(upHist, downHist);
                else
                {
                    delete upHist;
                    delete downHist;
                }
            }

            if (!channelShapeUncerts.empty())
                nonStackShapeUncerts.emplace(channel, std::move(channelShapeUncerts));
        }

        auto buildStandaloneChannelShapeUncerts = [&](const std::string &channel)
        {
            std::vector<std::pair<TH1D *, TH1D *>> channelShapeUncerts;
            auto existingIter = nonStackShapeUncerts.find(channel);
            if (existingIter != nonStackShapeUncerts.end())
            {
                int pairIndex = 0;
                for (const auto &[upHist, downHist] : existingIter->second)
                {
                    channelShapeUncerts.emplace_back(
                        cloneHistogramDetached(upHist, channel + "_standalone_up_" + std::to_string(pairIndex)),
                        cloneHistogramDetached(downHist, channel + "_standalone_down_" + std::to_string(pairIndex))
                    );
                    pairIndex++;
                }
            }

            auto isDataIter = isDataCopy.find(channel);
            if (isDataIter != isDataCopy.end() && isDataIter->second == 1)
                return channelShapeUncerts;

            const bool isSignalChannel = std::find(isSignal.begin(), isSignal.end(), channel) != isSignal.end();
            const bool isDDChannel = (isDDCopy[channel] == 1);
            const auto &campaignsToLoad = isSignalChannel ? signalCampaign : (isDDChannel ? ddCampaign : mcCampaign);
            const float scalingToUse = isSignalChannel ? signalScaling : (isDDChannel ? ddScaling : mcScaling);
            const std::string nominalWeightToUse = isSignalChannel ? MCWeight : (isDDChannel ? DDWeight : MCWeight);
            const auto &channelWeightSysts = isSignalChannel ? signalWeightSysts : (isDDChannel ? ddWeightSysts : weightSysts);
            const auto &channelShiftSystNames = isSignalChannel ? signalShiftSystNames : (isDDChannel ? ddShiftSystNames : shiftSystNames);

            std::vector<double> cropRange = {};
            if (needCrop[varName] == 1)
                cropRange = jsonConfig.at("cropedRange").at(varName).get<std::vector<double>>();

            if (existingIter == nonStackShapeUncerts.end())
            {
                for (const auto &[systName, variations] : channelWeightSysts)
                {
                    const std::string upWeightName = nominalWeightToUse + "_" + systName + "_up";
                    const std::string downWeightName = nominalWeightToUse + "_" + systName + "_down";
                    TH1D *upHist = loadChannelVariationTotal(inDir, campaignsToLoad, isDataCopy, processScale, needMerge, channel, varName, upWeightName, nominalWeightToUse, scalingToUse, varConfig, allowMissingAsZero, silentMissingSystForCrossEra, needCrop[varName] == 1, cropRange);
                    TH1D *downHist = loadChannelVariationTotal(inDir, campaignsToLoad, isDataCopy, processScale, needMerge, channel, varName, downWeightName, nominalWeightToUse, scalingToUse, varConfig, allowMissingAsZero, silentMissingSystForCrossEra, needCrop[varName] == 1, cropRange);
                    if (options.doNormalize)
                    {
                        if (upHist != nullptr && upHist->Integral() > 0.0)
                            upHist->Scale(1.0 / upHist->Integral());
                        if (downHist != nullptr && downHist->Integral() > 0.0)
                            downHist->Scale(1.0 / downHist->Integral());
                    }
                    if (upHist != nullptr && downHist != nullptr)
                        channelShapeUncerts.emplace_back(upHist, downHist);
                    else
                    {
                        delete upHist;
                        delete downHist;
                    }
                }
            }
            if (existingIter == nonStackShapeUncerts.end())
            {
                for (const auto &systName : channelShiftSystNames)
                {
                    const std::string upWeightName = nominalWeightToUse + "_" + systName + "_up";
                    const std::string downWeightName = nominalWeightToUse + "_" + systName + "_down";
                    TH1D *upHist = loadChannelVariationTotal(inDir, campaignsToLoad, isDataCopy, processScale, needMerge, channel, varName, upWeightName, nominalWeightToUse, scalingToUse, varConfig, allowMissingAsZero, silentMissingSystForCrossEra, needCrop[varName] == 1, cropRange);
                    TH1D *downHist = loadChannelVariationTotal(inDir, campaignsToLoad, isDataCopy, processScale, needMerge, channel, varName, downWeightName, nominalWeightToUse, scalingToUse, varConfig, allowMissingAsZero, silentMissingSystForCrossEra, needCrop[varName] == 1, cropRange);
                    if (options.doNormalize)
                    {
                        if (upHist != nullptr && upHist->Integral() > 0.0)
                            upHist->Scale(1.0 / upHist->Integral());
                        if (downHist != nullptr && downHist->Integral() > 0.0)
                            downHist->Scale(1.0 / downHist->Integral());
                    }
                    if (upHist != nullptr && downHist != nullptr)
                        channelShapeUncerts.emplace_back(upHist, downHist);
                    else
                    {
                        delete upHist;
                        delete downHist;
                    }
                }
            }

            if (isDDChannel)
            {
                auto ddVarIter = ddVariations.find(channel);
                auto nominalIter = histsNeeded.find(channel);
                if (ddVarIter != ddVariations.end() && nominalIter != histsNeeded.end() && nominalIter->second != nullptr)
                {
                    const double relVariation = ddVarIter->second;
                    if (relVariation > 0.0)
                    {
                        TH1D *upHist = cloneHistogramDetached(nominalIter->second, channel + "_DDvariation_up");
                        TH1D *downHist = cloneHistogramDetached(nominalIter->second, channel + "_DDvariation_down");
                        upHist->Scale(1.0 + relVariation);
                        downHist->Scale(std::max(0.0, 1.0 - relVariation));
                        if (options.doNormalize)
                        {
                            if (upHist->Integral() > 0.0)
                                upHist->Scale(1.0 / upHist->Integral());
                            if (downHist->Integral() > 0.0)
                                downHist->Scale(1.0 / downHist->Integral());
                        }
                        channelShapeUncerts.emplace_back(upHist, downHist);
                    }
                }
            }

            return channelShapeUncerts;
        };

        auto writeStandaloneSummaryObjects = [&](const std::string &rootFileName,
                                                 const std::string &channelKey,
                                                 TH1D *sourceHist,
                                                 const std::vector<std::pair<TH1D *, TH1D *>> &localShapeUncerts,
                                                 bool includeMCStat)
        {
            if (sourceHist == nullptr)
                return;

            std::unique_ptr<TFile> outFile(TFile::Open(rootFileName.c_str(), "UPDATE"));
            if (!outFile || outFile->IsZombie())
            {
                rdfWS_utility::messageERROR("plotHists", "Failed to open standalone summary ROOT file: " + rootFileName);
            }

            TH1D *standaloneHist = cloneHistogramDetached(sourceHist, channelKey);
            standaloneHist->SetDirectory(nullptr);
            outFile->cd();
            standaloneHist->Write(channelKey.c_str(), TObject::kOverwrite);

            TGraphAsymmErrors *uncertBand = buildStandaloneUncertaintyBand(standaloneHist, localShapeUncerts, includeMCStat);
            if (uncertBand != nullptr)
            {
                uncertBand->SetName((channelKey + "_uncert").c_str());
                uncertBand->Write(uncertBand->GetName(), TObject::kOverwrite);
            }

            delete standaloneHist;
            delete uncertBand;
        };

        // Nov 26, change the scaling in to a map
        // pHelper.drawStackHistWithRatio(histsNeeded, stackOrder, stackUp, stackDown, reOrder, ratioHists, options, mcScaling, colorScheme, channelLabels, drawHeader, drawText, {}, binLabels);
        pHelper.drawStackHistWithRatio(histsNeeded, stackOrder, stackUp, stackDown, reOrder, ratioHists, options, 1.0, colorScheme, channelLabels, drawHeader, drawText, {}, shapeUncerts, nonStackShapeUncerts, overlaySystConfig.includeMCStat, binLabels);

        const std::string standaloneRootFileName = plotName + "_objects.root";
        {
            std::unique_ptr<TFile> outFile(TFile::Open(standaloneRootFileName.c_str(), "RECREATE"));
            if (!outFile || outFile->IsZombie())
            {
                rdfWS_utility::messageERROR("plotHists", "Failed to create standalone summary ROOT file: " + standaloneRootFileName);
            }
        }

        TH1D *stackTotalHist = nullptr;
        for (const auto &stackChannel : stackOrder)
        {
            auto histIter = histsNeeded.find(stackChannel);
            if (histIter == histsNeeded.end() || histIter->second == nullptr)
                continue;
            if (stackTotalHist == nullptr)
            {
                stackTotalHist = cloneHistogramDetached(histIter->second, "stack_total_" + varName);
            }
            else
            {
                stackTotalHist->Add(histIter->second);
            }
        }
        if (stackTotalHist != nullptr)
        {
            writeStandaloneSummaryObjects(standaloneRootFileName, "stack_total", stackTotalHist, shapeUncerts, overlaySystConfig.includeMCStat);
            delete stackTotalHist;
        }

        for (const auto &[channel, hist] : histsNeeded)
        {
            const auto channelShapeUncerts = buildStandaloneChannelShapeUncerts(channel);
            writeStandaloneSummaryObjects(standaloneRootFileName, "channel_" + channel, hist, channelShapeUncerts, overlaySystConfig.includeMCStat);
            for (const auto &[upHist, downHist] : channelShapeUncerts)
            {
                delete upHist;
                delete downHist;
            }
        }

        for (auto &[upHist, downHist] : shapeUncerts)
        {
            delete upHist;
            delete downHist;
        }
        shapeUncerts.clear();
        for (auto &[channel, variations] : nonStackShapeUncerts)
        {
            for (auto &[upHist, downHist] : variations)
            {
                delete upHist;
                delete downHist;
            }
        }
        nonStackShapeUncerts.clear();
        for (auto &[histName, hist] : ratioHists)
        {
            delete hist;
        }
        ratioHists.clear();

        // delocate memory of copied hists
        for (auto &[histName, hist] : histsNeeded)
        {
            delete hist;
        }
        histsNeeded.clear();
    }

    return 0;
}



