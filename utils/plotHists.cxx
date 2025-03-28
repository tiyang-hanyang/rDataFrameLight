#include "PlotControl.h"
#include "HistControl.h"

#include "Utility.h"

#include <iostream>
#include <exception>
#include <string>
#include <map>
#include <vector>
#include <sstream>

#include "TH1D.h"

std::vector<std::map<std::string, double>> getStackUncert(std::vector<std::string> stackOrder, rdfWS_utility::JsonObject jsonConfig)
{
    if (stackOrder.size() == 0) return {};
    std::string systFile = jsonConfig.at("stackUncert");
    std::map<std::string, std::vector<double>> systJson = rdfWS_utility::readJson("plotHists", systFile);

    std::map<std::string, double> systUp, systDown;
    for (auto& [key, variation]: systJson)
    {
        systUp.emplace(key, variation[1]);
        systDown.emplace(key, variation[0]);
    }
    return {systUp, systDown};
}


int main(int argc, char *argv[])
{
    // read in the config
    if (argc < 2)
    {
        rdfWS_utility::messageERROR("plotHists", "No hist plot job json provided!");
    }
    std::string jsonPath = argv[1];
    // auto jsonConfig = rdfWS_utility::readJson("plotHists", jsonPath);
    rdfWS_utility::JsonObject jsonConfig(rdfWS_utility::readJson("plotHists", jsonPath), "Job Config");

    // parse basic info from json
    std::string jobName = jsonConfig.at("name");
    std::vector<std::string> runEra = jsonConfig.at("era");

    // formulate draw texts
    std::vector<std::string> drawText = jsonConfig.at("texts");
    std::string lumiPath = jsonConfig.at("lumiConfig");
    rdfWS_utility::JsonObject lumiConfig(rdfWS_utility::readJson("plotHists", lumiPath), "Lumi Config");
    float lumiValue = lumiConfig.at(runEra[0]);
    for (int i = 1; i < runEra.size(); i++)
    {
        float tempLumiValue = lumiConfig.at(runEra[i]);
        lumiValue += tempLumiValue;
    }
    for (int i = 0; i < drawText.size(); i++)
    {
        std::string text = drawText[i];
        if (text.find("%1.0f") != std::string::npos)
            drawText[i] = Form(text.c_str(), lumiValue);
    }

    // load channels
    std::vector<std::string> channels = jsonConfig.at("datasets");
    std::map<std::string, std::string> channelLabels = jsonConfig.at("datasetLabel");
    std::map<std::string, std::vector<std::string>> needMerge = jsonConfig.at("needMerge");
    std::vector<std::string> loadChannels;
    for (const auto &ch : channels)
    {
        if (needMerge.find(ch) == needMerge.end())
            loadChannels.push_back(ch);
        else
        {
            auto mergeList = needMerge[ch];
            for (auto mergeCh : mergeList)
            {
                loadChannels.push_back(mergeCh);
            }
        }
    }
    // for possible stack plots
    std::vector<std::string> stackOrder = jsonConfig.at("stackOrder");
    int reOrder = jsonConfig.at("reOrder");
    std::vector<std::string> numerator = jsonConfig.at("numerator");

    // treat styles
    const std::vector<std::string> isSignal = jsonConfig.at("isSignal");
    const std::map<std::string, int> isData = jsonConfig.at("isData");
    std::string dataWeight = jsonConfig.at("dataWeight");
    std::string MCWeight = jsonConfig.at("MCWeight");
    const std::map<std::string, int> colorScheme = jsonConfig.at("colorMapping");

    std::string inDir = jsonConfig.at("inDir");

    // for each variable, load the histograms and plot
    std::vector<std::string> variables = jsonConfig.at("varNames");
    std::map<std::string, int> needCrop = jsonConfig.at("needCrop");
    for (const auto &varName : variables)
    {
        // setup drawing options
        PlotContext options;
        std::vector<std::string> isSignalCopy = isSignal;
        options.isSignal = isSignalCopy;
        std::map<std::string, int> isDataCopy = isData;
        options.isData.push_back(isDataCopy);
        std::string varConfigPath = jsonConfig.at("varConfig");
        rdfWS_utility::JsonObject varJson(rdfWS_utility::readJson("plotHists", varConfigPath), "Var Config");
        auto varConfig = varJson.at(varName);
        options.xLabel = varConfig.at("label").get<std::string>();
        options.yLabel = {jsonConfig.at("yLabel"), jsonConfig.at("yRatioLabel")};
        options.doLog = jsonConfig.at("logPlot");
        options.xSize = jsonConfig.at("histXSize");
        options.ySize = jsonConfig.at("histYSize");

        // configure output
        std::string outputDir = jsonConfig.at("outDir");
        rdfWS_utility::creatingFolder("plotHists", outputDir);
        std::string plotName = outputDir + "/data_MC";
        for (auto era : runEra)
        {
            plotName += "_";
            plotName += era;
        }
        plotName += "_" + varName + "_" + jobName;
        if (needCrop[varName] == 1)
            plotName += "_crop";
        PlotControl pHelper(plotName);

        // read the histograms in
        HistControl histLoader;
        for (const auto &ch : loadChannels)
        {
            std::string histName = ch + "_" + varName + "_";
            if (isDataCopy[ch] == 1)
                histName += dataWeight;
            else
                histName += MCWeight;
            histLoader.loadHistogram(inDir + "_" + runEra[0] + "/" + histName + ".root", histName, ch, varName);
        }
        // add together when more than 1 eras
        for (int i = 1; i < runEra.size(); i++)
        {
            HistControl tempHistLoader;
            for (const auto &ch : loadChannels)
            {
                std::string histName = ch + "_" + varName + "_";
                if (isDataCopy[ch] == 1)
                    histName += dataWeight;
                else
                    histName += MCWeight;
                tempHistLoader.loadHistogram(inDir + "_" + runEra[i] + "/" + histName + ".root", histName, ch, varName);
            }
            histLoader = histLoader.addHistograms(tempHistLoader);
        }

        // take crops TODO
        if (needCrop[varName] == 1)
        {
            std::vector<int> cropRange = jsonConfig.at("cropedRange").at(varName);
            ;
            histLoader = histLoader.cropHistograms(cropRange[0], cropRange[1]);
        }

        // merge the plots according to the config
        for (const auto &[mergedCh, components] : needMerge)
        {
            histLoader.mergeHistograms(components, mergedCh);
        }

        // output for working plot
        auto histsNeeded = histLoader.getHists(channels);
        int doRatio = jsonConfig.at("doRatio");
        options.isData.push_back(std::map<std::string, int>{});
        for (auto key : numerator)
        {
            options.isData[1].emplace(key, isDataCopy[key]);
        }

        // allowing normalization, if no stack
        if (stackOrder.size() == 0)
        {
            if (jsonConfig.at("normalization").get<int>() == 1)
            {
                for (auto hist : histsNeeded)
                {
                    hist.second->Scale(1.0 / hist.second->Integral());
                }
            }
        }

        std::map<std::string, TH1D *> ratioHists = {};
        // currently only allowing ratio plot with respect to stack
        if (doRatio && stackOrder.size() > 0)
        {
            ratioHists = histLoader.getRatios(numerator, stackOrder);
        }

        // to add possible uncertainties
        auto stackUncert = getStackUncert(stackOrder, jsonConfig);
        std::map<std::string, double> stackUp = stackUncert[0];
        std::map<std::string, double> stackDown = stackUncert[1];
        pHelper.drawStackHistWithRatio(histsNeeded, stackOrder, stackUp, stackDown, reOrder, ratioHists, options, colorScheme, channelLabels, drawText);
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