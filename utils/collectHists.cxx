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

void prepareHist(rdfWS_utility::JsonObject jsonConfig, std::string variable, SampleControl samples)
{
    // extract dataset info
    std::vector<std::string> allChannels = jsonConfig.at("datasets");
    std::map<std::string, std::vector<std::string>> needMerge = jsonConfig.at("needMerge");
    std::map<std::string, int> isData = jsonConfig.at("isData");

    // extract var binning info
    std::string varConfigPath = jsonConfig.at("varConfig");
    rdfWS_utility::JsonObject varConfig(rdfWS_utility::readJson("collectHists", varConfigPath), "Var Config");
    // create binning object
    HistBinning *histBins = new HistBinning;
    histBins->nBins = std::stoi(varConfig.at(variable).at("nBins").get<std::string>());
    histBins->min = std::stof(varConfig.at(variable).at("min").get<std::string>());
    histBins->max = std::stof(varConfig.at(variable).at("max").get<std::string>());

    // extract cutflow
    CutControl histCut;
    std::vector<std::string> cutConfigList = jsonConfig.at("cuts");
    if (cutConfigList.size() > 0)
    {
        histCut = CutControl(cutConfigList[0]);
        for (int cutStep = 1; cutStep < cutConfigList.size(); cutStep++)
        {
            histCut = histCut + CutControl(cutConfigList[cutStep]);
        }
    }

    // readin lumi by era
    std::string lumiConfigPath = jsonConfig.at("lumiConfig");
    rdfWS_utility::JsonObject lumiConfig(rdfWS_utility::readJson("collectHists", lumiConfigPath), "Lumi Config");
    std::string era = jsonConfig.at("era");
    float lumiVal = lumiConfig.at(era);

    // out dir for storing histograms
    std::string outputDir = jsonConfig.at("outDir");
    // add the runEra infor together
    outputDir += "_";
    outputDir += era;
    outputDir += "/";

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
            std::vector<std::string> mcWeights = jsonConfig.at("MCweight");
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

        delete chDS;
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

    rdfWS_utility::JsonObject jsonConfig(rdfWS_utility::readJson("collectHists", jsonPath), "JO Config");
    // variables & samples
    std::vector<std::string> variables = jsonConfig.at("varNames");
    std::string sampleConfigPath = jsonConfig.at("sampleConfig");
    SampleControl samples(sampleConfigPath);
    for (const auto &varName : variables)
    {
        prepareHist(jsonConfig, varName, samples);
    }
    return 0;
}