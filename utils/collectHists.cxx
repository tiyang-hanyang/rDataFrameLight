#include "HistControl.h"
#include "SampleControl.h"
#include "CutControl.h"

#include "Utility.h"

#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <map>

#include "ROOT/RDataFrame.hxx"
#include "TChain.h"

void prepareHist(nlohmann::json jsonConfig, std::string variable, SampleControl samples)
{
    // extract dataset info
    std::vector<std::string> allChannels = jsonConfig["datasets"];
    std::map<std::string, std::vector<std::string>> needMerge = jsonConfig["needMerge"];
    std::map<std::string, int> isData = jsonConfig["isData"];

    // extract var binning info
    std::string varConfigPath = jsonConfig["varConfig"];
    auto varConfig = rdfWS_utility::readJson("collectHists", varConfigPath);
    // create binning object
    HistBinning *histBins = new HistBinning;
    histBins->nBins = std::stoi(std::string(varConfig[variable]["nBins"]));
    histBins->min = std::stof(std::string(varConfig[variable]["min"]));
    histBins->max = std::stof(std::string(varConfig[variable]["max"]));

    // extract cutflow
    CutControl histCut;
    std::vector<std::string> cutConfigList = jsonConfig["cuts"];
    if (cutConfigList.size() > 0)
    {
        histCut = CutControl(cutConfigList[0]);
        for (int cutStep = 1; cutStep < cutConfigList.size(); cutStep++)
        {
            histCut = histCut + CutControl(cutConfigList[cutStep]);
        }
    }

    // readin lumi by era
    std::string lumiConfigPath = jsonConfig["lumiConfig"];
    auto lumiConfig = rdfWS_utility::readJson("collectHists", lumiConfigPath);
    std::string era = jsonConfig["era"];
    float lumiVal = lumiConfig[era];

    // out dir for storing histograms
    std::string outputDir = jsonConfig["outDir"];

    // creating HistControl contrainer for all the channels
    HistControl varHistController;

    for (const auto &channel : allChannels)
    {
        // get file list
        std::vector<std::string> filePaths = {};
        if (needMerge.find(channel) != needMerge.end())
        {
            for (const auto &compName : needMerge[channel])
            {
                auto compPaths = samples.getFiles(compName);
                filePaths.insert(filePaths.end(), compPaths.begin(), compPaths.end());
            }
        }
        else
            filePaths = samples.getFiles(channel);

        if (filePaths.size() == 0)
        {
            rdfWS_utility::messageWARN("collectHists", "Sample list for " + channel + " is empty!");
            continue;
        }

        // create RooDataFrame and doing extraction
        TChain *chDS = new TChain("Events");
        for (auto filePath : filePaths)
        {
            chDS->Add(filePath.c_str());
        }
        ROOT::RDataFrame rdfDS(*chDS);
        ROOT::RDF::RNode rndDS(rdfDS);
        rndDS = histCut.applyCut(rndDS);

        // setup weight
        std::string weightName = "one";
        auto originalBrs = rndDS.GetColumnNames();
        if (std::find(originalBrs.begin(), originalBrs.end(), "one") == originalBrs.end())
            rndDS = rndDS.Define("one", "1");
        if (!isData[channel])
        {
            std::vector<std::string> mcWeights = jsonConfig["MCweight"];
            if (mcWeights.size() > 0)
            {
                weightName = "MCTotalWeight";
                if (std::find(originalBrs.begin(), originalBrs.end(), mcWeights[0]) == originalBrs.end())
                    rdfWS_utility::messageERROR("collectHists", "MC weight component " + mcWeights[0] + " not defined in channel " + channel + "! Please check your skim.");
                std::string MCWeightexp = "(" + mcWeights[0];
                for (int i = 1; i < mcWeights.size(); i++)
                {
                    std::string wName = mcWeights[i];
                    if (std::find(originalBrs.begin(), originalBrs.end(), wName) == originalBrs.end())
                        rdfWS_utility::messageERROR("collectHists", "MC weight component " + wName + " not defined in channel " + channel + "! Please check your skim.");
                    MCWeightexp += " * ";
                    MCWeightexp += wName;
                }
                MCWeightexp += " * " + std::to_string(lumiVal) + ")";
                rndDS = rndDS.Define("MCTotalWeight", MCWeightexp);
            }
        }

        // saving name would be datasetName+"_"+varName+"_"+weightName;
        varHistController.createHistogram(rndDS, channel, histBins, variable, weightName, outputDir);
    }

    delete histBins;
}

int main(int argc, char *argv[])
{
    ROOT::EnableImplicitMT();

    // readin the config
    if (argc < 2)
    {
        rdfWS_utility::messageERROR("collectHists", "No hist plot job json provided!");
    }
    std::string jsonPath = argv[1];

    auto jsonConfig = rdfWS_utility::readJson("collectHists", jsonPath);
    // variables & samples
    std::vector<std::string> variables = jsonConfig["varNames"];
    SampleControl samples(jsonConfig["sampleConfig"]);
    for (const auto &varName : variables)
    {
        prepareHist(jsonConfig, varName, samples);
    }
    return 0;
}