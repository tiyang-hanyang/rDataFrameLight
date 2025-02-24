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
    // auto jsonConfig = rdfWS_utility::readJson("plotHists", jsonPath);
    rdfWS_utility::JsonObject jsonConfig(rdfWS_utility::readJson("plotHists", jsonPath), "Job Config");

    // parse basic info from json
    std::string jobName = jsonConfig.at("name");
    std::string runEra = jsonConfig.at("era");

    // formulate draw texts
    std::vector<std::string> drawText = jsonConfig.at("texts");
    std::string lumiPath = jsonConfig.at("lumiConfig");
    rdfWS_utility::JsonObject lumiConfig(rdfWS_utility::readJson("plotHists", lumiPath), "Lumi Config");
    float lumiValue = lumiConfig.at(runEra);
    for (int i = 0; i < drawText.size(); i++)
    {
        std::string text = drawText[i];
        if (text.find("%1.0f") != std::string::npos)
            drawText[i] = Form(text.c_str(), lumiValue);
    }

    // load channels
    std::vector<std::string> channels = jsonConfig.at("datasets");
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

    std::map<std::string, int> isData = jsonConfig.at("isData");
    std::string dataWeight = jsonConfig.at("dataWeight");
    std::string MCWeight = jsonConfig.at("MCWeight");

    std::string inDir = jsonConfig.at("inDir");

    // for each variable, load the histograms and plot
    std::vector<std::string> variables = jsonConfig.at("varNames");
    std::map<std::string, int> needCrop = jsonConfig.at("needCrop");
    for (const auto &varName : variables)
    {
        // setup drawing options
        PlotContext options;
        options.isData.push_back(isData);
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
            // TEST
            // std::vector<int> cropRange = jsonConfig["cropedRange"][varName];
            // std::cout << "flag 1" << std::endl;
            // auto cropRange = rdfWS_utility::parseJson<std::vector<int>>(jsonConfig["cropedRange"], "job option", varName, "plotHists");
            // std::cout << "flag 2" << std::endl;
            // auto cropRange2 = jsonConfig.at("cropedRange").at(varName).get<std::vector<int>>();
            // std::cout << "flag 3" << std::endl;
            // auto cropRange3 = jsonConfig.at("cropedRange").at("aaaa").get<std::vector<int>>();
            // std::cout << "flag 4" << std::endl;
            // auto cropRangeWrong = rdfWS_utility::parseJson<std::vector<int>>(jsonConfig["cropedRange"], "job option", "aaaa", "plotHists");
            // std::cout << "flag 5" << std::endl;
            std::vector<int> cropRange = jsonConfig.at("cropedRange").at(varName);;
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
                if (jsonConfig.at("normalization").get<int>() == 1)
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