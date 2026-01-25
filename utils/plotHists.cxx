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

// function to read systmatic variables
// if non provided, will return empty vector and plot nominal-only
std::vector<std::map<std::string, double>> getStackUncert(std::vector<std::string> stackOrder, rdfWS_utility::JsonObject jsonConfig)
{
    if (stackOrder.size() == 0)
        return {};
    std::string systFile = jsonConfig.at("systConfig");
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
    if (configFile.contains("jobType"))
    {
        std::string jobType = configFile.at("jobType").get<std::string>();
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

    // Nov 26
    // add the data-driven contribution, these contribution is from data but used similar to the MC
    int withDD(0);
    std::vector<std::string> ddCampaign;
    std::vector<std::string> ddChannels;
    if (jsonConfig.contains("dd_era"))
    {
        withDD=1;
        ddCampaign = jsonConfig.at("dd_era").get<std::vector<std::string>>();
        if (jsonConfig.contains("dd_datasets"))
            ddChannels = jsonConfig.at("dd_datasets").get<std::vector<std::string>>();
        else
        {
            rdfWS_utility::messageWARN("plotHists", "DD era added, but dd_datasets not assigned, using nonprompt as default");
            ddChannels.emplace_back("nonprompt");
        }
    }

    // get lumi-value of the plot
    std::string lumiPath = jsonConfig.at("lumiConfig");
    rdfWS_utility::JsonObject lumiConfig(rdfWS_utility::readJson("plotHists", lumiPath), "Lumi Config");
    float lumiValue = lumiConfig.at(runEra[0]);
    for (int i = 1; i < runEra.size(); i++)
    {
        float tempLumiValue = lumiConfig.at(runEra[i]);
        lumiValue += tempLumiValue;
    }

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
    float mclumi = lumiConfig.at(mcCampaign[0]);
    for (int i = 1; i < mcCampaign.size(); i++)
    {
        float tempLumiValue = lumiConfig.at(mcCampaign[i]);
        mclumi += tempLumiValue;
    }
    auto mcScaling = lumiValue / mclumi;
    std::cout << "mc lumi scaling value:" << mcScaling << std::endl;

    // DD rescaling
    float ddScaling(0.0);
    if (withDD)
    {
        float ddlumi = lumiConfig.at(ddCampaign[0]);
        for (int i = 1; i < ddCampaign.size(); i++)
        {
            float tempLumiValue = lumiConfig.at(ddCampaign[i]);
            ddlumi += tempLumiValue;
        }
        ddScaling = lumiValue / ddlumi;
        std::cout << "data-driven lumi scaling value:" << ddScaling << std::endl;
    }

    // load channels
    std::vector<std::string> channels = jsonConfig.at("datasets");
    std::map<std::string, std::string> channelLabels = jsonConfig.at("datasetLabel");
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
        options.doLog = jsonConfig.at("doLogPlot");
        options.doNormalize = jsonConfig.at("doNormPlot");

        std::vector<std::string> isSignalCopy = isSignal;
        options.isSignal = isSignalCopy;
        std::map<std::string, int> isDataCopy = isData;
        options.isData.push_back(isDataCopy);
        std::string varConfigPath = jsonConfig.at("varConfig");
        rdfWS_utility::JsonObject varJson(rdfWS_utility::readJson("plotHists", varConfigPath), "Var Config");
        auto varConfig = varJson.at(varName);
        options.xLabel = varConfig.at("label").get<std::string>();

        // add default values for y-axis labels
        std::string yUpperLabel = "Events";
        if (options.doNormalize) yUpperLabel = "a.u.";
        std::string yLowerLabel = "Data / MC";
        if (configFile.contains("yLabel"))
        {
            yUpperLabel = configFile.at("yLabel").get<std::string>();
        }
        if (configFile.contains("yRatioLabe"))
        {
            yLowerLabel = configFile.at("yRatioLabel").get<std::string>();
        }
        options.yLabel = {yUpperLabel, yLowerLabel};
        options.xSize = jsonConfig.at("histXSize");
        options.ySize = jsonConfig.at("histYSize");

        std::vector<std::string> binLabels = varConfig.at("binLabels");
        // need to truncate if have cropping
        if (needCrop[varName] == 1 && binLabels.size()>0)
        {
            std::cout << "print debugging" << std::endl;
            std::cout << "orig bin label size: " << binLabels.size() << std::endl;
            for (auto& labelText: binLabels) std::cout << "\t" << labelText << std::endl;
            std::vector<double> cropRange = jsonConfig.at("cropedRange").at(varName);
            std::string strOrigMax = varConfig.at("max");
            double origMax = std::stod(strOrigMax);
            std::cout << "original max: " << origMax << std::endl;
            std::string strOrigMin = varConfig.at("min");
            double origMin = std::stod(strOrigMin);
            std::cout << "original min: " << origMin << std::endl;
            std::string strOrigBins = varConfig.at("nBins");
            int origBins = std::stoi(strOrigBins);
            std::cout << "original bins: " << origBins << std::endl;
            double originalWidth = (origMax - origMin) / origBins;
            std::cout << "bin width: " << originalWidth << std::endl;
            int first_bin = static_cast<int>((cropRange[0] - origMin) / originalWidth);
            std::cout << "new bin start from: " << first_bin << std::endl;
            int last_bin = static_cast<int>((cropRange[1] - origMin) / originalWidth);
            std::cout << "new bin up to: " << last_bin << std::endl;
            if (first_bin < 0) first_bin = 0;
            if (last_bin >= origBins) last_bin = origBins - 1;
            binLabels = std::vector<std::string>(binLabels.begin() + first_bin, binLabels.begin() + last_bin + 1);
            std::cout << "new bin label size: " << binLabels.size() << std::endl;
            for (auto& labelText: binLabels) std::cout << "\t" << labelText << std::endl;
            std::cout << "truncate bin labels done!" << std::endl;
        }

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
            histLoader.loadHistogram(inDir + "_" + runEra[0] + "/" + histName + ".root", histName, ch, 1.0, varName);
        }

        // first mc era: mcCampaign[0]
        HistControl mcHistLoader;
        int withMC=0;
        for (const auto &ch : loadChannels)
        {
            std::string histName = ch + "_" + varName + "_";
            if (isDataCopy[ch] == 1)
                continue;
            else
                histName += MCWeight;
            mcHistLoader.loadHistogram(inDir + "_" + mcCampaign[0] + "/" + histName + ".root", histName, ch, mcScaling, varName);
            withMC++;
        }
        if (withMC>0)
            histLoader = histLoader.addHistograms(mcHistLoader);

        // Nov 26, add DD contribution
        if (withDD)
        {
            HistControl DDHistLoader;
            // dd does not merge histograms, so only need to directly read from the config
            for (const auto &ch : ddChannels)
            {
                // now just hard coding this variable, later maybe I should also try to contain the MC backgrounds
                std::string histName = ch + "_" + varName + "_" + MCWeight;
                DDHistLoader.loadHistogram(inDir + "_" + ddCampaign[0] + "/" + histName + ".root", histName, ch, ddScaling, varName);
            }
            histLoader = histLoader.addHistograms(DDHistLoader);

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
                else
                    continue;
                tempHistLoader.loadHistogram(inDir + "_" + runEra[i] + "/" + histName + ".root", histName, ch, 1.0, varName);
                histNum++;
            }
            if (histNum>0) 
                histLoader = histLoader.addHistograms(tempHistLoader);
        }

        // mc campaigns
        for (int i = 1; i < mcCampaign.size(); i++)
        {
            HistControl tempHistLoader;
            for (const auto &ch : loadChannels)
            {
                std::string histName = ch + "_" + varName + "_";
                if (isDataCopy[ch] == 1)
                    continue;
                else
                    histName += MCWeight;
                tempHistLoader.loadHistogram(inDir + "_" + mcCampaign[i] + "/" + histName + ".root", histName, ch, mcScaling, varName);
            }
            histLoader = histLoader.addHistograms(tempHistLoader);
        }

        // Nov 26, add DD contribution
        if (withDD)
        {
            for (int i = 1; i < mcCampaign.size(); i++)
            {
                HistControl tempDDHistLoader;
                // dd does not merge histograms, so only need to directly read from the config
                for (const auto &ch : ddChannels)
                {
                    // now just hard coding this variable, later maybe I should also try to contain the MC backgrounds
                    std::string histName = ch + "_" + varName + "_" + MCWeight;
                    tempDDHistLoader.loadHistogram(inDir + "_" + ddCampaign[i] + "/" + histName + ".root", histName, ch, ddScaling, varName);
                }
                histLoader = histLoader.addHistograms(tempDDHistLoader);
            }
        }

        // take crops TODO
        if (needCrop[varName] == 1)
        {
            std::vector<double> cropRange = jsonConfig.at("cropedRange").at(varName);
            histLoader = histLoader.cropHistograms(cropRange[0], cropRange[1]);
        }

        // merge the plots according to the config
        for (const auto &[mergedCh, components] : needMerge)
        {
            histLoader.mergeHistograms(components, mergedCh);
        }

        // output for working plot
        // Nov 26, dd channel must be added.
        std::vector<std::string> allChannels = channels;
        if (withDD)
        {
            allChannels.insert(allChannels.end(), ddChannels.begin(), ddChannels.end());
        }
        auto histsNeeded = histLoader.getHists(allChannels);
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
        auto stackUncert = getStackUncert(stackOrder, jsonConfig);
        std::map<std::string, double> stackUp, stackDown;
        if (stackUncert.size() == 2)
        {
            stackUp = stackUncert[0];
            stackDown = stackUncert[1];
        }
        // Nov 26, change the scaling in to a map
        // pHelper.drawStackHistWithRatio(histsNeeded, stackOrder, stackUp, stackDown, reOrder, ratioHists, options, mcScaling, colorScheme, channelLabels, drawHeader, drawText, {}, binLabels);
        pHelper.drawStackHistWithRatio(histsNeeded, stackOrder, stackUp, stackDown, reOrder, ratioHists, options, 1.0, colorScheme, channelLabels, drawHeader, drawText, {}, binLabels);
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
