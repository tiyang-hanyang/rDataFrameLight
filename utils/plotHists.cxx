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

int main(int argc, char *argv[])
{
    // read in the config
    if (argc < 2)
    {
        rdfWS_utility::messageERROR("plotHists", "No hist plot job json provided!");
    }
    std::string jsonPath = argv[1];
    auto jsonConfig = rdfWS_utility::readJson("plotHists", jsonPath);

    // parse basic info from json
    std::string jobName = jsonConfig["name"];
    std::string runEra = jsonConfig["era"];

    // formulate draw texts
    std::vector<std::string> drawText = jsonConfig["texts"];
    float lumiValue = rdfWS_utility::readJson("plotHists", jsonConfig["lumiConfig"])[runEra];
    for (int i = 0; i < drawText.size(); i++)
    {
        std::string text = drawText[i];
        if (text.find("%1.0f") != std::string::npos)
            drawText[i] = Form(text.c_str(), lumiValue);
    }

    // load channels
    std::vector<std::string> channels = jsonConfig["datasets"];
    std::map<std::string, std::vector<std::string>> needMerge = jsonConfig["needMerge"];
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
    std::vector<std::string> stackOrder = jsonConfig["stackOrder"];
    int reOrder = jsonConfig["reOrder"];
    std::vector<std::string> numerator = jsonConfig["numerator"];

    std::map<std::string, int> isData = jsonConfig["isData"];
    std::string dataWeight = jsonConfig["dataWeight"];
    std::string MCWeight = jsonConfig["MCWeight"];

    std::string inDir = jsonConfig["inDir"];

    // for each variable, load the histograms and plot
    std::vector<std::string> variables = jsonConfig["varNames"];
    std::map<std::string, int> needCrop = jsonConfig["needCrop"];
    for (const auto &varName : variables)
    {
        // setup drawing options
        PlotContext options;
        options.isData.push_back(isData);
        std::string varConfigPath = jsonConfig["varConfig"];
        auto varConfig = rdfWS_utility::readJson("plotHists", varConfigPath)[varName];
        options.xLabel = varConfig["label"];
        options.yLabel = {jsonConfig["yLabel"], jsonConfig["yRatioLabel"]};
        options.doLog = jsonConfig["logPlot"];
        options.xSize = jsonConfig["histXSize"];
        options.ySize = jsonConfig["histYSize"];

        // configure output
        std::string outputDir = jsonConfig["outDir"];
        rdfWS_utility::creatingFolder("plotHists", outputDir);
        std::string plotName = outputDir + "/data_MC_" + runEra + "_" + varName + "_" + jobName;
        if (needCrop[varName] == 1)
            plotName += "_crop";
        PlotControl pHelper(plotName);

        // read the histograms in
        HistControl histLoader;
        for (const auto &ch : loadChannels)
        {
            std::string histName = ch + "_" + varName + "_";
            if (isData[ch] == 1)
                histName += dataWeight;
            else
                histName += MCWeight;
            histLoader.loadHistogram(inDir + "/" + histName + ".root", histName, ch, varName);
        }

        // take crops TODO
        if (needCrop[varName] == 1)
        {   
            std::vector<int> cropRange = jsonConfig["cropedRange"][varName];
            histLoader = histLoader.cropHistograms(cropRange[0], cropRange[1]);
        }

        // merge the plots according to the config
        for (const auto &[mergedCh, components] : needMerge)
        {
            histLoader.mergeHistograms(components, mergedCh);
        }

        // output for working plot
        auto histsNeeded = histLoader.getHists(channels);
        int doRatio = jsonConfig["doRatio"];
        if (doRatio)
        {
            options.isData.push_back(std::map<std::string, int>{});
            for (auto key : numerator)
            {
                options.isData[1].emplace(key, isData[key]);
            }
            auto ratioHists = histLoader.getRatios(numerator, stackOrder);
            pHelper.drawStackHistWithRatio(histsNeeded, stackOrder, reOrder, ratioHists, options, drawText);
            for (auto &[histName, hist] : ratioHists)
            {
                delete hist;
            }
            ratioHists.clear();
        }
        else
        {
            if (stackOrder.size() > 0)
                pHelper.drawStackHist(histsNeeded, stackOrder, reOrder, options, drawText);
            else
            {
                if (jsonConfig["normalization"] == 1)
                {
                    for (auto hist : histsNeeded)
                    {
                        hist.second->Scale(1.0 / hist.second->Integral());
                    }
                }

                pHelper.drawHist(histsNeeded, options, drawText);
            }
        }

        // delocate memory of copied hists
        for (auto &[histName, hist] : histsNeeded)
        {
            delete hist;
        }
        histsNeeded.clear();
    }

    return 0;
}