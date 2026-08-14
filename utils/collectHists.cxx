#include "HistControl.h"
#include "SampleControl.h"
#include "CutControl.h"

#include "Utility.h"

#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <utility>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <limits>
#include <memory>

#include "ROOT/RDataFrame.hxx"
#include "TChain.h"
#include "TH2D.h"

// processed result of the sample list of one channel.
// @validFilePaths: samples that contain valid Events.
// @totalWeight: the summed genWeightSum include even samples without valid Events inside.
struct ChannelSample
{
    std::vector<std::string> validFilePaths;
    double totalWeight = 0.0;
};

struct ShiftSystematic
{
    std::string name;
    std::vector<std::string> targets;
    std::string nominalWeightBranch;
    std::string upWeightBranch;
    std::string downWeightBranch;
    std::string upWeightExpression;
    std::string downWeightExpression;
};

struct WeightSystematicSpec
{
    std::string name;
    std::string nominalWeight;
    std::string upWeight;
    std::string downWeight;
};

// variable histogram structure from json
// enables multiple dimensions
struct VarHistInfo
{
    int dimension = 1;
    std::string xVar;
    std::string yVar;
    std::string zVar;
    HistBinning xBin;
    HistBinning yBin;
    HistBinning zBin;

    std::string histVarName() const
    {
        if (dimension == 1)
            return xVar;
        if (dimension == 2)
            return xVar + "_vs_" + yVar;
        if (dimension == 3)
            return xVar + "_vs_" + yVar + "_vs_" + zVar;
        return xVar;
    }

    VarHistInfo resolve_shift_syst(const std::function<std::string(const std::string &)> &resolver) const
    {
        VarHistInfo shiftedVarInfo = *this;
        shiftedVarInfo.xVar = resolver(xVar);
        shiftedVarInfo.yVar = resolver(yVar);
        shiftedVarInfo.zVar = resolver(zVar);
        return shiftedVarInfo;
    }
};

struct HistRequest
{
    std::string name;
    std::string weightVar;
    VarHistInfo varInfo;
};

struct BookedHistHandle
{
    virtual ~BookedHistHandle() = default;
    virtual TH1 *cloneHist(const std::string &histName) = 0;
};

template <typename HistT>
struct TypedBookedHistHandle : BookedHistHandle
{
    ROOT::RDF::RResultPtr<HistT> hist;

    explicit TypedBookedHistHandle(ROOT::RDF::RResultPtr<HistT> bookedHist)
        : hist(std::move(bookedHist))
    {
    }

    TH1 *cloneHist(const std::string &histName) override
    {
        TH1 *histClone = static_cast<TH1 *>(hist.GetPtr()->Clone(histName.c_str()));
        histClone->SetDirectory(0);
        return histClone;
    }
};

struct BookedHist
{
    std::string name;
    std::unique_ptr<BookedHistHandle> handle;
};

namespace
{
    std::map<std::string, std::vector<std::string>> loadSystAliases()
    {
        rdfWS_utility::JsonObject aliasConfig(
            rdfWS_utility::readJson("collectHists", "json/general_config/syst_NP.json"),
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
}

double getAdditionalScale(rdfWS_utility::JsonObject jsonConfig, const std::string &channel)
{
    double scaleFactor = 1.0;
    if (!jsonConfig.contains("scale"))
        return scaleFactor;

    std::vector<std::string> scaleConfigList = jsonConfig.at("scale");
    for (const auto &scaleConfigPath : scaleConfigList)
    {
        rdfWS_utility::JsonObject scaleConfig(rdfWS_utility::readJson("collectHists", scaleConfigPath), "Scale Config");
        if (scaleConfig.contains(channel))
        {
            double channelScale = scaleConfig.at(channel);
            scaleFactor *= channelScale;
        }
    }
    return scaleFactor;
}

ChannelSample get_valid_file_list(
    SampleControl &samples,
    const std::string &channel,
    bool requireGenWeightSum)
{
    ChannelSample sample;
    auto filePaths = samples.getFiles(channel);
    if (filePaths.empty())
    {
        rdfWS_utility::messageWARN("collectHists", "Channel " + channel + " has no files listed in SampleControl. Skip histogram creation.");
        return sample;
    }

    for (const auto &tempF : filePaths)
    {
        std::unique_ptr<TFile> f_temp(TFile::Open(tempF.c_str(), "READ"));
        if (!f_temp || f_temp->IsZombie())
        {
            rdfWS_utility::messageWARN("collectHists", "Cannot open file " + tempF);
            continue;
        }

        if (f_temp->GetListOfKeys()->FindObject("genWeightSum"))
        {
            auto hsum = dynamic_cast<TH1D *>(f_temp->Get("genWeightSum"));
            if (!hsum)
            {
                rdfWS_utility::messageWARN("collectHists", "no genWeightSum hist in " + tempF);
                if (requireGenWeightSum)
                    continue;
            }
            else
            {
                sample.totalWeight += hsum->Integral();
            }
        }
        else if (requireGenWeightSum)
        {
            rdfWS_utility::messageWARN("collectHists", "no genWeightSum hist in " + tempF + ", skip.");
            continue;
        }

        if (f_temp->GetListOfKeys()->FindObject("Events"))
            sample.validFilePaths.push_back(tempF);
    }

    return sample;
}


// July 6
double get_MC_scale(
    rdfWS_utility::JsonObject &lumiConfig,
    rdfWS_utility::JsonObject &XSConfig,
    const std::string &era,
    const std::string &channel,
    double totalWeight,
    const std::vector<std::string> &mcWeights)
{
    double scaleFactor = 1.0;
    const double lumiVal = lumiConfig.at(era);

    if (std::find(mcWeights.begin(), mcWeights.end(), "weight_XS") != mcWeights.end())
    {
        scaleFactor = lumiVal;
    }
    else
    {
        const double XSval = XSConfig.at(channel);
        if (totalWeight > 0.0)
            scaleFactor = lumiVal * XSval * 1000.0 / totalWeight;
        else
            rdfWS_utility::messageWARN("collectHists", "Non-positive total MC weight for " + channel + ". Keep MC scale as 1.");
    }
    return scaleFactor;
}

ROOT::RDF::TH1DModel makeTH1DModelFromVarInfo(const VarHistInfo &varInfo)
{
    if (varInfo.dimension != 1)
        throw std::runtime_error("TH1D model build requires a 1D variable.");

    const HistBinning *binning = &varInfo.xBin;
    const std::string &title = varInfo.xVar;
    if (!binning)
        throw std::runtime_error("Cannot build histogram model from null binning.");

    if (!binning->varBins.empty())
    {
        if (binning->nBins != -1 && static_cast<int>(binning->varBins.size()) != binning->nBins + 1)
            throw std::runtime_error("Variable binning size does not match nBins + 1.");
        return ROOT::RDF::TH1DModel("templateHistModel", title.c_str(), static_cast<int>(binning->varBins.size()) - 1, binning->varBins.data());
    }

    if (binning->nBins == -1)
        throw std::runtime_error("Lazy histogram booking requires explicit binning.");

    return ROOT::RDF::TH1DModel("templateHistModel", title.c_str(), binning->nBins, binning->min, binning->max);
}

ROOT::RDF::TH2DModel makeTH2DModelFromVarInfo(const VarHistInfo &varInfo)
{
    if (varInfo.dimension != 2)
        throw std::runtime_error("TH2D model build requires a 2D variable.");

    const HistBinning &xBin = varInfo.xBin;
    const HistBinning &yBin = varInfo.yBin;
    const std::string title = varInfo.histVarName();
    const bool xVarBins = !xBin.varBins.empty();
    const bool yVarBins = !yBin.varBins.empty();

    const auto makeFixedEdges = [](const HistBinning &bin, const std::string &axisName)
    {
        if (bin.nBins <= 0 || bin.max <= bin.min)
            throw std::runtime_error("Fixed " + axisName + " binning requires positive nBins and max > min.");
        std::vector<double> edges(static_cast<std::size_t>(bin.nBins) + 1);
        const double width = (bin.max - bin.min) / static_cast<double>(bin.nBins);
        for (int index = 0; index <= bin.nBins; ++index)
            edges[static_cast<std::size_t>(index)] = bin.min + index * width;
        return edges;
    };

    if (xVarBins && yVarBins)
    {
        return ROOT::RDF::TH2DModel(
            "templateHistModel",
            title.c_str(),
            static_cast<int>(xBin.varBins.size()) - 1,
            xBin.varBins.data(),
            static_cast<int>(yBin.varBins.size()) - 1,
            yBin.varBins.data()
        );
    }
    if (xVarBins && !yVarBins)
    {
        const auto yEdges = makeFixedEdges(yBin, "y-axis");
        return ROOT::RDF::TH2DModel(
            "templateHistModel",
            title.c_str(),
            static_cast<int>(xBin.varBins.size()) - 1,
            xBin.varBins.data(),
            static_cast<int>(yEdges.size()) - 1,
            yEdges.data()
        );
    }

    if (!xVarBins && yVarBins)
    {
        const auto xEdges = makeFixedEdges(xBin, "x-axis");
        return ROOT::RDF::TH2DModel(
            "templateHistModel",
            title.c_str(),
            static_cast<int>(xEdges.size()) - 1,
            xEdges.data(),
            static_cast<int>(yBin.varBins.size()) - 1,
            yBin.varBins.data()
        );
    }

    return ROOT::RDF::TH2DModel(
        "templateHistModel",
        title.c_str(),
        xBin.nBins,
        xBin.min,
        xBin.max,
        yBin.nBins,
        yBin.min,
        yBin.max
    );
}

BookedHist book_hist(ROOT::RDF::RNode node, const HistRequest &request)
{
    if (request.varInfo.dimension == 1)
    {
        auto histModel = makeTH1DModelFromVarInfo(request.varInfo);
        return {
            request.name,
            std::make_unique<TypedBookedHistHandle<TH1D>>(
                node.Histo1D(histModel, request.varInfo.xVar, request.weightVar)
            )
        };
    }

    if (request.varInfo.dimension == 2)
    {
        auto histModel = makeTH2DModelFromVarInfo(request.varInfo);
        return {
            request.name,
            std::make_unique<TypedBookedHistHandle<TH2D>>(
                node.Histo2D(histModel, request.varInfo.xVar, request.varInfo.yVar, request.weightVar)
            )
        };
    }

    throw std::runtime_error("Unsupported histogram dimension: " + std::to_string(request.varInfo.dimension));
}

void fold_overflow_underflow_1d(TH1 *hist)
{
    const int nBinsX = hist->GetNbinsX();
    std::vector<double> contents(nBinsX + 1, 0.0);
    std::vector<double> errors2(nBinsX + 1, 0.0);

    for (int ix = 0; ix <= nBinsX + 1; ++ix)
    {
        const int targetX = std::min(std::max(ix, 1), nBinsX);
        const double content = hist->GetBinContent(ix);
        const double error = hist->GetBinError(ix);
        contents[targetX] += content;
        errors2[targetX] += error * error;
    }

    for (int ix = 0; ix <= nBinsX + 1; ++ix)
    {
        hist->SetBinContent(ix, 0.0);
        hist->SetBinError(ix, 0.0);
    }

    for (int ix = 1; ix <= nBinsX; ++ix)
    {
        hist->SetBinContent(ix, contents[ix]);
        hist->SetBinError(ix, std::sqrt(errors2[ix]));
    }
}

void fold_overflow_underflow_2d(TH1 *hist)
{
    const int nBinsX = hist->GetNbinsX();
    const int nBinsY = hist->GetNbinsY();
    std::vector<std::vector<double>> contents(nBinsX + 1, std::vector<double>(nBinsY + 1, 0.0));
    std::vector<std::vector<double>> errors2(nBinsX + 1, std::vector<double>(nBinsY + 1, 0.0));

    for (int ix = 0; ix <= nBinsX + 1; ++ix)
    {
        for (int iy = 0; iy <= nBinsY + 1; ++iy)
        {
            const int targetX = std::min(std::max(ix, 1), nBinsX);
            const int targetY = std::min(std::max(iy, 1), nBinsY);
            const double content = hist->GetBinContent(ix, iy);
            const double error = hist->GetBinError(ix, iy);
            contents[targetX][targetY] += content;
            errors2[targetX][targetY] += error * error;
        }
    }

    for (int ix = 0; ix <= nBinsX + 1; ++ix)
    {
        for (int iy = 0; iy <= nBinsY + 1; ++iy)
        {
            hist->SetBinContent(ix, iy, 0.0);
            hist->SetBinError(ix, iy, 0.0);
        }
    }

    for (int ix = 1; ix <= nBinsX; ++ix)
    {
        for (int iy = 1; iy <= nBinsY; ++iy)
        {
            hist->SetBinContent(ix, iy, contents[ix][iy]);
            hist->SetBinError(ix, iy, std::sqrt(errors2[ix][iy]));
        }
    }
}

void fold_overflow_underflow(TH1 *hist)
{
    if (hist->GetDimension() == 1)
    {
        fold_overflow_underflow_1d(hist);
        return;
    }
    if (hist->GetDimension() == 2)
    {
        fold_overflow_underflow_2d(hist);
        return;
    }
}

void save_channel_histograms(
    const std::vector<TH1 *> &hists,
    const std::string &histName,
    const std::string &outputDir)
{
    if (hists.empty())
        return;

    if (!std::filesystem::exists(outputDir))
        std::filesystem::create_directories(outputDir);

    TFile *f_save = TFile::Open((outputDir + "/" + histName + ".root").c_str(), "RECREATE");
    if (!f_save || f_save->IsZombie())
    {
        delete f_save;
        throw std::runtime_error("Failed to create output ROOT file: " + outputDir + "/" + histName + ".root");
    }
    TH1 *saveHist = static_cast<TH1 *>(hists[0]->Clone(histName.c_str()));
    saveHist->SetDirectory(0);
    for (auto iter = hists.begin() + 1; iter != hists.end(); ++iter)
        saveHist->Add(*iter);
    fold_overflow_underflow(saveHist);

    saveHist->Write("");
    f_save->Close();
    delete f_save;
    delete saveHist;
}

void materialize_and_save(
    std::vector<BookedHist> &bookedHists,
    double scaleFactor,
    const std::string &outputDir)
{
    for (auto &booked : bookedHists)
    {
        TH1 *histClone = booked.handle->cloneHist(booked.name);
        histClone->Scale(scaleFactor);
        save_channel_histograms({histClone}, booked.name, outputDir);
        delete histClone;
    }
}

std::string buildWeightExpressionWithReplacement(
    const std::vector<std::string> &weights,
    const std::string &nominalWeight,
    const std::string &replacementWeight
)
{
    if (weights.empty())
        return "";

    std::string expression = "(";
    for (size_t i = 0; i < weights.size(); ++i)
    {
        if (i > 0)
            expression += " * ";
        expression += (weights[i] == nominalWeight) ? replacementWeight : weights[i];
    }
    expression += ")";
    return expression;
}

std::string buildShiftedWeightBranchName(
    const std::string &nominalWeight,
    const std::string &systName,
    bool isUp)
{
    if (nominalWeight == systName)
        return nominalWeight + (isUp ? "_up" : "_down");

    std::string systSuffix = systName;
    if (nominalWeight == "JVMweight" && systSuffix.rfind("CMS_", 0) == 0)
        systSuffix = systSuffix.substr(4);
    if (nominalWeight == "btag_weight" && systSuffix.rfind("btag_", 0) == 0)
        systSuffix = systSuffix.substr(5);
    return nominalWeight + "_" + systSuffix + (isUp ? "_up" : "_down");
}

std::vector<std::string> get_mc_weights(
    rdfWS_utility::JsonObject jsonConfig,
    const std::string &nominalWeightName)
{
    if (nominalWeightName == "one")
        return {};

    const std::string weightKey = (nominalWeightName == "DDTotalWeight") ? "DDweight" : "MCweight";
    if (!jsonConfig.contains(weightKey))
        rdfWS_utility::messageERROR("collectHists", weightKey + " is missing in the config.");

    std::vector<std::string> baseWeights = jsonConfig.at(weightKey).get<std::vector<std::string>>();
    if (baseWeights.empty())
        rdfWS_utility::messageERROR("collectHists", weightKey + " is empty in the config.");

    return baseWeights;
}

std::vector<WeightSystematicSpec> buildWeightSystematics(
    rdfWS_utility::JsonObject jsonConfig,
    const std::string &nominalWeightName)
{
    std::vector<WeightSystematicSpec> specs;
    const auto baseWeights = get_mc_weights(jsonConfig, nominalWeightName);
    if (!jsonConfig.contains("weightSyst"))
        return specs;

    const auto aliases = loadSystAliases();
    auto rawConfig = jsonConfig.at("weightSyst").get<nlohmann::json>();
    if (!rawConfig.is_object())
        rdfWS_utility::messageERROR("collectHists", "weightSyst must be a map.");

    for (auto it = rawConfig.begin(); it != rawConfig.end(); ++it)
    {
        const std::string rawName = it.key();
        if (!it.value().is_string())
            rdfWS_utility::messageERROR("collectHists", "weightSyst entry " + rawName + " must be a string nominal weight.");

        const std::string nominalWeight = it.value().get<std::string>();
        bool hasNominalWeight = false;
        for (const auto &weight : baseWeights)
        {
            if (weight == nominalWeight)
            {
                hasNominalWeight = true;
                break;
            }
        }
        if (!hasNominalWeight)
        {
            rdfWS_utility::messageWARN("collectHists", "weightSyst " + rawName + " nominal weight " + nominalWeight + " is not in nominal weight list; skip this weight systematic.");
            continue;
        }

        for (const auto &systName : expandAliasOrName(rawName, aliases))
        {
            specs.push_back({
                systName,
                nominalWeight,
                buildShiftedWeightBranchName(nominalWeight, systName, true),
                buildShiftedWeightBranchName(nominalWeight, systName, false)
            });
        }
    }
    return specs;
}

std::vector<ShiftSystematic> buildShiftSystematics(
    rdfWS_utility::JsonObject jsonConfig,
    const std::string &nominalWeightName
)
{
    std::vector<ShiftSystematic> shiftSysts;
    const auto baseWeights = get_mc_weights(jsonConfig, nominalWeightName);
    const auto aliases = loadSystAliases();

    std::map<std::string, std::vector<std::string>> shiftConfig;
    if (jsonConfig.contains("shiftSyst"))
    {
        auto rawConfig = jsonConfig.at("shiftSyst").get<nlohmann::json>();
        if (!rawConfig.is_object())
            rdfWS_utility::messageERROR("collectHists", "shiftSyst must be a map.");
        for (auto it = rawConfig.begin(); it != rawConfig.end(); ++it)
        {
            const std::string rawName = it.key();
            if (!it.value().is_array())
                rdfWS_utility::messageERROR("collectHists", "shiftSyst entry " + rawName + " must be a list of target branches.");
            const auto targets = it.value().get<std::vector<std::string>>();
            for (const auto &name : expandAliasOrName(rawName, aliases))
            {
                if (shiftConfig.find(name) != shiftConfig.end())
                    rdfWS_utility::messageERROR("collectHists", "Duplicate shiftSyst after alias expansion: " + name);
                shiftConfig.emplace(name, targets);
            }
        }
    }
    std::map<std::string, std::vector<std::string>> shiftWeightConfig;
    if (jsonConfig.contains("shiftWeightSyst"))
    {
        auto rawConfig = jsonConfig.at("shiftWeightSyst").get<nlohmann::json>();
        if (!rawConfig.is_object())
            rdfWS_utility::messageERROR("collectHists", "shiftWeightSyst must be a map.");
        for (auto it = rawConfig.begin(); it != rawConfig.end(); ++it)
        {
            const std::string rawName = it.key();
            if (!it.value().is_array())
                rdfWS_utility::messageERROR("collectHists", "shiftWeightSyst entry " + rawName + " must be a list.");
            const auto config = it.value().get<std::vector<std::string>>();
            const auto expandedNames = expandAliasOrName(rawName, aliases);
            if (expandedNames.size() > 1 && config.size() != 1)
            {
                rdfWS_utility::messageERROR("collectHists", "shiftWeightSyst alias " + rawName + " only supports shorthand form [nominalWeight].");
            }
            for (const auto &name : expandedNames)
            {
                if (shiftWeightConfig.find(name) != shiftWeightConfig.end())
                    rdfWS_utility::messageERROR("collectHists", "Duplicate shiftWeightSyst after alias expansion: " + name);
                if (config.size() == 1)
                {
                    shiftWeightConfig.emplace(name, std::vector<std::string>{
                        config[0],
                        buildShiftedWeightBranchName(config[0], name, true),
                        buildShiftedWeightBranchName(config[0], name, false)
                    });
                }
                else
                {
                    shiftWeightConfig.emplace(name, config);
                }
            }
        }
    }

    std::set<std::string> systNames;
    for (const auto &[systName, targets] : shiftConfig)
        systNames.insert(systName);

    for (const auto &systName : systNames)
    {
        ShiftSystematic syst;
        syst.name = systName;
        if (shiftConfig.find(systName) != shiftConfig.end())
            syst.targets = shiftConfig.at(systName);
        if (syst.targets.empty())
        {
            rdfWS_utility::messageERROR("collectHists", "shiftSyst " + syst.name + " has no target branch.");
            exit(1);
        }
        syst.upWeightExpression = buildWeightExpressionWithReplacement(baseWeights, "", "");
        syst.downWeightExpression = syst.upWeightExpression;
        auto weightIt = shiftWeightConfig.find(syst.name);
        if (weightIt != shiftWeightConfig.end())
        {
            const auto &weights = weightIt->second;
            if (weights.size() != 3)
            {
                rdfWS_utility::messageERROR("collectHists", "shiftWeightSyst " + syst.name + " must be [nominalWeight, upWeight, downWeight].");
                exit(1);
            }
            bool foundNominal = false;
            for (const auto &weight : baseWeights)
            {
                if (weight == weights[0])
                {
                    foundNominal = true;
                    break;
                }
            }
            if (!foundNominal)
            {
                rdfWS_utility::messageERROR("collectHists", "shiftWeightSyst " + syst.name + " nominal weight " + weights[0] + " is not in nominal weight list.");
                exit(1);
            }
            syst.nominalWeightBranch = weights[0];
            syst.upWeightBranch = weights[1];
            syst.downWeightBranch = weights[2];
            syst.upWeightExpression = buildWeightExpressionWithReplacement(baseWeights, weights[0], weights[1]);
            syst.downWeightExpression = buildWeightExpressionWithReplacement(baseWeights, weights[0], weights[2]);
        }
        shiftSysts.push_back(syst);
    }
    return shiftSysts;
}

std::string build_shifted_target_name(
    const std::string &target,
    const std::string &systName,
    bool isUp)
{
    const std::string direction = isUp ? "up" : "down";
    const auto pos = target.find('_');
    if (pos == std::string::npos)
        return target + "_" + systName + "_" + direction;
    return target.substr(0, pos) + "_" + systName + "_" + direction + "_" + target.substr(pos + 1);
}

std::string build_shifted_target_expression(
    const std::string &target,
    const std::string &systName,
    bool isUp)
{
    if (systName == "JER_corr")
        return target + "/Jet_JER_corr*Jet_JER_corr_" + std::string(isUp ? "up" : "down");

    return target + (isUp ? "*(1+" : "*(1-") + systName + ")";
}

// July 6
// test the validity of the systematics
bool check_weight_systematic_validity(
    ROOT::RDF::RNode rndDS,
    const WeightSystematicSpec &syst,
    const std::vector<std::string> &baseWeights,
    std::string &reason)
{
    if (std::find(baseWeights.begin(), baseWeights.end(), syst.nominalWeight) == baseWeights.end())
    {
        reason = "nominal weight " + syst.nominalWeight + " is not used in nominal weight list";
        return false;
    }

    if (!rndDS.HasColumn(syst.nominalWeight))
    {
        reason = "nominal weight column " + syst.nominalWeight + " is missing";
        return false;
    }

    if (!rndDS.HasColumn(syst.upWeight) || !rndDS.HasColumn(syst.downWeight))
    {
        reason = "systematic weight columns " + syst.upWeight + "/" + syst.downWeight + " are missing";
        return false;
    }
    reason.clear();
    return true;
}

// July 6
bool check_shift_systematic_validity(
    ROOT::RDF::RNode rndDS,
    const ShiftSystematic &syst,
    bool isUp,
    std::string &reason
)
{
    for (const auto &target : syst.targets)
    {
        if (!rndDS.HasColumn(target))
        {
            reason = "target column " + target + " is missing";
            return false;
        }

        if (syst.name == "JER_corr")
        {
            if (!rndDS.HasColumn("Jet_JER_corr"))
            {
                reason = "required column Jet_JER_corr is missing";
                return false;
            }
            const std::string jerShiftColumn = "Jet_JER_corr_" + std::string(isUp ? "up" : "down");
            if (!rndDS.HasColumn(jerShiftColumn))
            {
                reason = "required column " + jerShiftColumn + " is missing";
                return false;
            }
        }
        else if (!rndDS.HasColumn(syst.name))
        {
            reason = "systematic column " + syst.name + " is missing";
            return false;
        }
    }

    if (!syst.nominalWeightBranch.empty())
    {
        if (!rndDS.HasColumn(syst.nominalWeightBranch))
        {
            reason = "nominal weight column " + syst.nominalWeightBranch + " is missing";
            return false;
        }
        if (!rndDS.HasColumn(syst.upWeightBranch) || !rndDS.HasColumn(syst.downWeightBranch))
        {
            reason = "systematic weight columns " + syst.upWeightBranch + "/" + syst.downWeightBranch + " are missing";
            return false;
        }
    }

    reason.clear();
    return true;
}

//////////
// July 6
// getting hist bins, must exist
std::unique_ptr<HistBinning> loadHistBinning(
    rdfWS_utility::JsonObject &varConfig,
    const std::string &variable
)
{
    auto varEntry = varConfig.at(variable);

    auto histBins = std::make_unique<HistBinning>();

    const bool hasNbins = varEntry.contains("nBins");
    const bool hasVarBinning = varEntry.contains("binning");
    const bool hasFixedMinMax = varEntry.contains("min") && varEntry.contains("max");
    if (!hasNbins || (!hasVarBinning && !hasFixedMinMax))
        rdfWS_utility::messageERROR("collectHists", "Variable " + variable + " must define either [nBins + binning] or [nBins + min + max].");

    // readin nBins
    histBins->nBins = std::stoi(varEntry.at("nBins").get<std::string>());
    if (histBins->nBins <= 0)
        rdfWS_utility::messageERROR("collectHists", "Variable " + variable + " must define a positive nBins.");

    // readin variable binning
    if (hasVarBinning)
    {
        histBins->varBins = varEntry.at("binning").get<std::vector<double>>();
        if (static_cast<int>(histBins->varBins.size()) != histBins->nBins + 1)
            rdfWS_utility::messageERROR("collectHists", "Variable binning size for " + variable + " must equal nBins + 1.");
        for (size_t i = 1; i < histBins->varBins.size(); ++i)
        {
            if (histBins->varBins[i] <= histBins->varBins[i - 1])
                rdfWS_utility::messageERROR("collectHists", "Variable binning edges for " + variable + " must be strictly increasing.");
        }
        histBins->min = histBins->varBins.front();
        histBins->max = histBins->varBins.back();
    }
    // readin fixed binning
    else
    {
        histBins->min = std::stof(varEntry.at("min").get<std::string>());
        histBins->max = std::stof(varEntry.at("max").get<std::string>());
        if (!std::isfinite(histBins->min) || !std::isfinite(histBins->max))
            rdfWS_utility::messageERROR("collectHists", "Fixed binning min/max for " + variable + " must be finite.");
        if (histBins->max <= histBins->min)
            rdfWS_utility::messageERROR("collectHists", "Fixed binning for " + variable + " requires max > min.");
    }
    return histBins;
}

std::vector<VarHistInfo> buildVarHistInfos(rdfWS_utility::JsonObject jsonConfig)
{
    std::vector<VarHistInfo> varInfos;
    const std::string varConfigPath = jsonConfig.at("varConfig");
    rdfWS_utility::JsonObject varConfig(rdfWS_utility::readJson("collectHists", varConfigPath), "Var Config");

    if (jsonConfig.contains("varNames"))
    {
        for (const auto &variable : jsonConfig.at("varNames").get<std::vector<std::string>>())
        {
            auto histBins = loadHistBinning(varConfig, variable);
            VarHistInfo varInfo;
            varInfo.dimension = 1;
            varInfo.xVar = variable;
            varInfo.xBin = *histBins;
            varInfos.push_back(varInfo);
        }
    }

    if (jsonConfig.contains("varNames_2D"))
    {
        auto rawVarNames2D = jsonConfig.at("varNames_2D").get<nlohmann::json>();
        if (!rawVarNames2D.is_array())
            rdfWS_utility::messageERROR("collectHists", "varNames_2D must be a list of [xVar, yVar] pairs.");

        for (const auto &entry : rawVarNames2D)
        {
            if (!entry.is_array() || entry.size() != 2)
                rdfWS_utility::messageERROR("collectHists", "Each varNames_2D entry must be [xVar, yVar].");

            const std::string xVar = entry.at(0).get<std::string>();
            const std::string yVar = entry.at(1).get<std::string>();
            auto xBins = loadHistBinning(varConfig, xVar);
            auto yBins = loadHistBinning(varConfig, yVar);

            VarHistInfo varInfo;
            varInfo.dimension = 2;
            varInfo.xVar = xVar;
            varInfo.yVar = yVar;
            varInfo.xBin = *xBins;
            varInfo.yBin = *yBins;
            varInfos.push_back(varInfo);
        }
    }

    if (varInfos.empty())
        rdfWS_utility::messageERROR("collectHists", "At least one variable must be provided in varNames or varNames_2D.");

    return varInfos;
}

// collect job type:
// 0 = MC by default
// 1 = DD only when "isDD": 1 (or true) is explicitly provided
// 2 = Data only when "isData": 1 (or true) is explicitly provided
int get_collect_job_type(rdfWS_utility::JsonObject jsonConfig)
{
    bool isDD = false;
    if (jsonConfig.contains("isDD"))
    {
        auto rawValue = jsonConfig.at("isDD").get<nlohmann::json>();
        if (rawValue.is_boolean())
            isDD = rawValue.get<bool>();
        else if (rawValue.is_number_integer())
            isDD = rawValue.get<int>() == 1;
        else
            rdfWS_utility::messageERROR("collectHists", "isDD must be 1/0 or true/false.");
    }

    bool isData = false;
    if (jsonConfig.contains("isData"))
    {
        auto rawValue = jsonConfig.at("isData").get<nlohmann::json>();
        if (rawValue.is_boolean())
            isData = rawValue.get<bool>();
        else if (rawValue.is_number_integer())
            isData = rawValue.get<int>() == 1;
        else
            rdfWS_utility::messageERROR("collectHists", "isData must be 1/0 or true/false.");
    }

    if (isDD && isData)
        rdfWS_utility::messageERROR("collectHists", "collectHists expects a single dataset type per job; do not enable both isData and isDD.");
    if (isDD)
        return 1;
    if (isData)
        return 2;
    return 0;
}

std::string build_nominal_weight_expression(
    rdfWS_utility::JsonObject jsonConfig,
    const std::string &nominalWeightName)
{
    if (nominalWeightName == "one")
        return "1.0";

    const std::string weightKey = (nominalWeightName == "MCTotalWeight") ? "MCweight" : "DDweight";
    if (!jsonConfig.contains(weightKey))
        rdfWS_utility::messageERROR("collectHists", weightKey + " is missing in the config.");

    std::vector<std::string> weights = jsonConfig.at(weightKey).get<std::vector<std::string>>();
    if (weights.empty())
        rdfWS_utility::messageERROR("collectHists", weightKey + " is empty in the config.");

    std::string expression = "(" + weights[0];
    for (int i = 1; i < static_cast<int>(weights.size()); i++)
    {
        expression += " * ";
        expression += weights[i];
    }
    expression += ")";
    return expression;
}

// final channel list for collect:
// - keeps the top-level dataset order
// - applies --channels filtering on top-level names
// - expands merged entries into their component channels
// second entry is the extra output subdir; empty means standalone
std::vector<std::pair<std::string, std::string>> get_collect_channels(
    rdfWS_utility::JsonObject jsonConfig,
    const std::set<std::string> &selectedChannels)
{
    std::vector<std::pair<std::string, std::string>> channels;
    std::vector<std::string> datasets = jsonConfig.at("datasets");
    std::set<std::string> topLevelChannels(datasets.begin(), datasets.end());

    for (const auto &selected : selectedChannels)
    {
        if (topLevelChannels.find(selected) == topLevelChannels.end())
        {
            rdfWS_utility::messageERROR("collectHists", "Selected channel " + selected + " is not a top-level dataset.");
        }
    }

    std::map<std::string, std::vector<std::string>> needMerge;
    if (jsonConfig.contains("needMerge"))
    {
        auto rawConfig = jsonConfig.at("needMerge").get<nlohmann::json>();
        if (!rawConfig.is_object())
            rdfWS_utility::messageERROR("collectHists", "needMerge must be a map from merged channel to component channel list.");
        for (auto it = rawConfig.begin(); it != rawConfig.end(); ++it)
        {
            if (!it.value().is_array())
                rdfWS_utility::messageERROR("collectHists", "needMerge entry " + it.key() + " must be a list.");
            needMerge[it.key()] = it.value().get<std::vector<std::string>>();
        }
    }

    for (const auto &dataset : datasets)
    {
        if (!selectedChannels.empty() && selectedChannels.find(dataset) == selectedChannels.end())
            continue;

        auto mergeIt = needMerge.find(dataset);
        if (mergeIt == needMerge.end())
        {
            channels.push_back({dataset, ""});
            continue;
        }

        for (const auto &component : mergeIt->second)
            channels.push_back({component, dataset});
    }

    return channels;
}

void process_channel(
    const std::string &histName,
    const std::vector<std::string> &filesWithEvents,
    CutControl &histCut,
    const VarHistInfo &varInfo,
    double scaleFactor,
    const std::vector<std::string> &baseWeights,
    const std::string &nominalWeightName,
    const std::string &nominalWeightExpr,
    const std::string &outputDir,
    const std::vector<WeightSystematicSpec> &weightSystematics,
    const std::vector<ShiftSystematic> &shiftSystematics
)
{
    // get into the file interaction
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame rdfDS("Events", filesWithEvents);

    // build the nominal node
    ROOT::RDF::RNode nomRndDS(rdfDS);
    std::set<std::string> preDefinedStepNames;
    nomRndDS = histCut.applyDefineOnly(
        nomRndDS,
        nullptr,
        [&](const std::string &name)
        {
            preDefinedStepNames.insert(name);
        });
    nomRndDS = histCut.applyCutSkippingSteps(nomRndDS, preDefinedStepNames);
    nomRndDS = nomRndDS.Define(nominalWeightName, nominalWeightExpr);

    // nominal and weight systematics share the same nominal-cut branch
    std::vector<BookedHist> bookedHists;
    bookedHists.push_back(book_hist(nomRndDS, {
        histName + "_" + nominalWeightName,
        nominalWeightName,
        varInfo
    }));

    // weight syst booking
    for (const auto &weightSyst : weightSystematics)
    {
        std::string skipReason;
        const bool canApplyWeight = check_weight_systematic_validity(nomRndDS, weightSyst, baseWeights, skipReason);
        if (!canApplyWeight)
        {
            rdfWS_utility::messageWARN("collectHists", "Weight systematic " + weightSyst.name + " is not available for " + histName + " (" + skipReason + "); skip this variation.");
            continue;
        }

        for (const auto &direction : std::vector<std::pair<std::string, std::string>>{
                 {"up", weightSyst.upWeight},
                 {"down", weightSyst.downWeight}})
        {
            const std::string systWeightName = nominalWeightName + "_" + weightSyst.name + "_" + direction.first;
            nomRndDS = nomRndDS.Define(
                systWeightName,
                buildWeightExpressionWithReplacement(baseWeights, weightSyst.nominalWeight, direction.second)
            );
            bookedHists.push_back(book_hist(nomRndDS, {
                histName + "_" + systWeightName,
                systWeightName,
                varInfo
            }));
        }
    }

    materialize_and_save(bookedHists, scaleFactor, outputDir);

    // non-weight syst hists booking (including redefinition on cuts)
    std::vector<BookedHist> bookedShiftHists;
    ROOT::RDF::RNode shiftRndDS(rdfDS);
    for (const auto &shiftSyst : shiftSystematics)
    {
        for (const auto &direction : std::vector<std::pair<std::string, bool>>{{"up", true}, {"down", false}})
        {
            const std::string shiftWeightName = nominalWeightName + "_" + shiftSyst.name + "_" + direction.first;
            const std::string suffix = "__" + shiftSyst.name + "_" + direction.first;
            const std::string passName = "pass_all" + suffix;
            std::string skipReason;
            const bool canApplyShift = check_shift_systematic_validity(shiftRndDS, shiftSyst, direction.second, skipReason);
            if (!canApplyShift)
            {
                rdfWS_utility::messageWARN(
                    "collectHists",
                    "Systematic " + shiftSyst.name + " is not available for " + histName
                    + " (" + skipReason + "); skip "
                    + direction.first + " variation."
                );
                continue;
            }

            std::map<std::string, std::string> baseReplacements;
            for (const auto &target : shiftSyst.targets)
            {
                const std::string shiftedName = build_shifted_target_name(target, shiftSyst.name, direction.second);
                const std::string shiftedExpr = build_shifted_target_expression(target, shiftSyst.name, direction.second);
                shiftRndDS = shiftRndDS.Define(shiftedName, shiftedExpr);
                baseReplacements[target] = shiftedName;
            }

            shiftRndDS = histCut.applySuffixedCutTags(shiftRndDS, suffix, baseReplacements, passName);
            shiftRndDS = shiftRndDS.Define(
                shiftWeightName,
                "(" + (direction.second ? shiftSyst.upWeightExpression : shiftSyst.downWeightExpression) + ")*(" + passName + ")"
            );

            VarHistInfo shiftedVarInfo = varInfo.resolve_shift_syst(
                [&](const std::string &varName)
                {
                    return histCut.resolveSuffixedColumnName(varName, suffix, baseReplacements);
                });
            bookedShiftHists.push_back(book_hist(shiftRndDS, {
                histName + "_" + shiftWeightName,
                shiftWeightName,
                shiftedVarInfo
            }));
        }
    }

    materialize_and_save(bookedShiftHists, scaleFactor, outputDir);
    ROOT::DisableImplicitMT();
}



void prepareHist(
    rdfWS_utility::JsonObject jsonConfig,
    const VarHistInfo &varInfo,
    SampleControl samples,
    const std::set<std::string> &selectedChannels
)
{
    const int jobType = get_collect_job_type(jsonConfig);
    const bool jobIsMC = (jobType == 0);
    const auto channels = get_collect_channels(jsonConfig, selectedChannels);

    // extract cutflow
    CutControl histCut;
    std::vector<std::string> cutConfigList = jsonConfig.at("cutConfig");
    if (cutConfigList.size() > 0)
    {
        histCut = CutControl(cutConfigList[0]);
        for (int cutStep = 1; cutStep < cutConfigList.size(); cutStep++)
        {
            histCut = histCut + CutControl(cutConfigList[cutStep]);
        }
    }

    std::string lumiConfigPath = jsonConfig.at("lumiConfig");
    rdfWS_utility::JsonObject lumiConfig(rdfWS_utility::readJson("collectHists", lumiConfigPath), "Lumi Config");
    std::string era = jsonConfig.at("era");

    std::string XSConfigPath = "json/XS/Run3.json";
    if (jsonConfig.contains("XSConfig"))
    {
        std::string XSloadPath = jsonConfig.at(std::string("XSConfig"));
        XSConfigPath = XSloadPath;
    }
    if (XSConfigPath == "")
        XSConfigPath = "json/XS/Run3.json";
    rdfWS_utility::JsonObject XSConfig(rdfWS_utility::readJson("collectHists", XSConfigPath), "XS Config");

    // out dir for storing histograms
    std::string outputDir = jsonConfig.at("outDir");
    // add the runEra infor together
    outputDir += "_" + era + "/";

    std::string nominalWeightName = "one";
    if (jobIsMC)
        nominalWeightName = "MCTotalWeight";
    else if (jobType == 1)
        nominalWeightName = "DDTotalWeight";
    std::string nominalWeightExpr = build_nominal_weight_expression(jsonConfig, nominalWeightName);

    std::vector<std::string> baseWeights;
    auto weightSystematics = std::vector<WeightSystematicSpec>{};
    auto shiftSystematics = std::vector<ShiftSystematic>{};
    if (nominalWeightName != "one")
    {
        baseWeights = get_mc_weights(jsonConfig, nominalWeightName);
        weightSystematics = buildWeightSystematics(jsonConfig, nominalWeightName);
        shiftSystematics = buildShiftSystematics(jsonConfig, nominalWeightName);
    }

    // processing by channel
    for (const auto &[channel, outputSubdir] : channels)
    {
        auto channelSample = get_valid_file_list(samples, channel, jobIsMC);
        if (channelSample.validFilePaths.empty())
        {
            rdfWS_utility::messageWARN("collectHists", "Channel " + channel + " has no valid ROOT files with Events tree after file existence checks. Skip histogram creation.");
            continue;
        }

        double scaleFactor = 1.0;
        if (jobIsMC)
            scaleFactor = get_MC_scale(lumiConfig, XSConfig, era, channel, channelSample.totalWeight, baseWeights);
        scaleFactor *= getAdditionalScale(jsonConfig, channel);

        rdfWS_utility::messageINFO("collectHists", "Starting channel " + channel);
        const std::string histName = channel + "_" + varInfo.histVarName();
        const std::string channelOutputDir = outputDir + (outputSubdir.empty() ? "" : outputSubdir + "/");
        process_channel(
            histName,
            channelSample.validFilePaths,
            histCut,
            varInfo,
            scaleFactor,
            baseWeights,
            nominalWeightName,
            nominalWeightExpr,
            channelOutputDir,
            weightSystematics,
            shiftSystematics
        );
    }
}

int main(int argc, char *argv[])
{
    // input menu
    // collectHists <json path> [--channels <channel1> <channel2> ...]
    if (argc < 2)
    {
        rdfWS_utility::messageERROR("collectHists", "No hist plot job json provided!");
    }
    std::string jsonPath = argv[1];
    std::set<std::string> selectedChannels;
    for (int i = 2; i < argc; i++)
    {
        std::string arg = argv[i];
        if (arg == "--channels")
            continue;
        selectedChannels.insert(arg);
    }

    rdfWS_utility::JsonObject jsonConfig(rdfWS_utility::readJson("collectHists", jsonPath), "JO Config");
    // job type check, better to have to avoid confusion
    if (jsonConfig.contains("jobType"))
    {
        std::string jobType = jsonConfig.at("jobType").get<std::string>();
        if (jobType != "collectingHists") 
        {
            rdfWS_utility::messageERROR("collectHists", "The jobType of your config is not collectingHists! Please check again. Running ceases...");
            exit(1);
        }
    }

    // variables & samples
    std::string sampleConfigPath = jsonConfig.at("sampleConfig");
    SampleControl samples(sampleConfigPath);
    auto varInfos = buildVarHistInfos(jsonConfig);
    for (const auto &varInfo : varInfos)
    {
        prepareHist(jsonConfig, varInfo, samples, selectedChannels);
    }
    return 0;
}
