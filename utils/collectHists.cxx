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

// Jan 19
// function to enable evaluating the lumi scaling at the plotting
// Then we could have the flexible input number of samples, else we need to determine before skimming
double gettingTotalWeight(const std::vector<std::string> &fileLists)
{
    double totalWeight(0.0);
    for (auto tempF : fileLists)
    {
        std::unique_ptr<TFile> f_temp(TFile::Open(tempF.c_str(), "READ"));
        if (!f_temp || f_temp->IsZombie())
        {
            rdfWS_utility::messageWARN("collectHists", "Cannot open file " + tempF);
            continue;
        }
        auto hsum = dynamic_cast<TH1D *>(f_temp->Get("genWeightSum"));
        if (!hsum)
        {
            rdfWS_utility::messageWARN("collectHists", "no genWeightSum hist in " + tempF);
            continue;
        }
        totalWeight += hsum->Integral();
    }
    return totalWeight;
}

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
    std::vector<std::string> cutConfigList = jsonConfig.at("cutConfig");
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

    // and also XS
    nlohmann::json XS_values;
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
    outputDir += "_";
    outputDir += era;
    outputDir += "/";

    // creating HistControl contrainer for all the channels
    HistControl varHistController;

    // common MC weight
    std::vector<std::string> mcWeights = jsonConfig.at("MCweight");
    std::string MCWeightexp = "";
    if (mcWeights.size() > 0)
    {
        MCWeightexp = "(" + mcWeights[0];
        for (int i = 1; i < mcWeights.size(); i++)
        {
            std::string wName = mcWeights[i];
            MCWeightexp += " * ";
            MCWeightexp += wName;
        }
        MCWeightexp += ")";
    }

    for (const auto &channel : allChannels)
    {
        // in case of merging histogram
        if (needMerge.find(channel) != needMerge.end())
        {

            std::vector<TH1D *> hists;
            // if any of the component is MC, then should be all MC, else should be all data
            std::string weightName = "one";
            // do one by one, need to compute the
            ROOT::EnableImplicitMT();
            for (const auto &compName : needMerge[channel])
            {

                // one of channel inside the merging scheme
                auto compPaths = samples.getFiles(compName);

                // scale
                double scaleFactor = 1.0;
                // should still looking for the whole channel's status, as if the sum is MC, then component is also MC
                if (!isData[channel])
                {
                    // Jan 20
                    // add backward compatibility, if we already have the weight_XS instead of the genWeight, we could try with only scaled by lumiVal instead of the following
                    // in this case, we have the weight_XS = genWeight * XS * 1000 / sum(genWeight) already from the beginning
                    if (MCWeightexp.find("weight_XS") != std::string::npos)
                    {
                        scaleFactor = lumiVal;
                    }
                    else
                    {
                        // need to get XSval
                        double XSval = XSConfig.at(compName);
                        // total weights
                        double totalWeight = gettingTotalWeight(compPaths);
                        // the get lumi-XS scaling
                        if (totalWeight > 0.0)
                            scaleFactor = lumiVal * XSval * 1000.0 / totalWeight;
                    }
                }

                std::vector<std::string> hasEventList;
                for (auto tempF : compPaths)
                {
                    std::unique_ptr<TFile> f_temp(TFile::Open(tempF.c_str(), "READ"));
                    if (!f_temp || f_temp->IsZombie())
                    {
                        rdfWS_utility::messageWARN("collectHists", "Cannot open file " + tempF);
                        continue;
                    }
                    if (f_temp->GetListOfKeys()->FindObject("Events"))
                    {
                        hasEventList.push_back(tempF);
                    }
                }

                if (hasEventList.size() == 0)
                    continue;

                // getting rdataframe for extracting histograms
                ROOT::RDataFrame rdfDS("Events", hasEventList);
                ROOT::RDF::RNode rndDS(rdfDS);
                rndDS = histCut.applyCut(rndDS);
                rndDS = rndDS.Define("one", "1");
                if (!isData[channel])
                {
                    rndDS = rndDS.Define("MCTotalWeight", MCWeightexp);
                    weightName = "MCTotalWeight";
                }

                // create histogram without saving, then using additional method to save
                // WARNING: the return object of creatHistogram is different from the ones read from the histogram TFiles, the container inside the HistController own the pointer as well, do not delete the pointer manually, or segfault will be triggered.
                // We are still considering if there is anything we could do to solve this
                auto componentHist = varHistController.createHistogram(rndDS, compName, histBins, variable, weightName, outputDir, 0);
                componentHist->Scale(scaleFactor);
                hists.push_back(componentHist);
            }

            if (hists.size() == 0)
            {
                rdfWS_utility::messageWARN("collectHists", channel+" has no events survive the selection, skip...");
                continue;
            }

            // summing histogram, save, and release memory (only needed in case of merging)
            std::string histName = channel + "_" + variable + "_" + weightName;
            ROOT::DisableImplicitMT();
            if (!std::filesystem::exists(outputDir))
                std::filesystem::create_directories(outputDir);
            TFile *f_save = new TFile((outputDir + "/" + histName + ".root").c_str(), "RECREATE");
            TH1D *saveHist = (TH1D *)hists[0]->Clone(histName.c_str());

            std::cout << saveHist->Integral() << std::endl;
            for (auto iter = hists.begin() + 1; iter < hists.end(); iter++)
            {
                saveHist->Add(*iter);
            }

            saveHist->Write("");

            f_save->Close();
        }
        else
        {
            auto filePaths = samples.getFiles(channel);

            double scaleFactor = 1.0;
            if (!isData[channel])
            {
                // Jan 20
                // add backward compatibility, if we already have the weight_XS instead of the genWeight, we could try with only scaled by lumiVal instead of the following
                if (MCWeightexp.find("weight_XS") != std::string::npos)
                {
                    scaleFactor = lumiVal;
                }
                else
                {
                    double XSval = XSConfig.at(channel);
                    double totalWeight = gettingTotalWeight(filePaths);
                    if (totalWeight > 0.0)
                        scaleFactor = lumiVal * XSval * 1000 / totalWeight;
                }
            }

            std::vector<std::string> hasEventList;
            for (auto tempF : filePaths)
            {
                std::unique_ptr<TFile> f_temp(TFile::Open(tempF.c_str(), "READ"));
                if (!f_temp || f_temp->IsZombie())
                {
                    rdfWS_utility::messageWARN("collectHists", "Cannot open file " + tempF);
                    continue;
                }
                if (f_temp->GetListOfKeys()->FindObject("Events"))
                {
                    hasEventList.push_back(tempF);
                }
            }

            std::string weightName = "one";
            ROOT::RDataFrame rdfDS("Events", hasEventList);
            ROOT::RDF::RNode rndDS(rdfDS);
            rndDS = histCut.applyCut(rndDS);
            rndDS = rndDS.Define("one", "1");
            if (!isData[channel])
            {
                rndDS = rndDS.Define("MCTotalWeight", MCWeightexp + "*" + std::to_string(scaleFactor));
                weightName = "MCTotalWeight";
            }

            varHistController.createHistogram(rndDS, channel, histBins, variable, weightName, outputDir);
        }
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
    std::vector<std::string> variables = jsonConfig.at("varNames");
    std::string sampleConfigPath = jsonConfig.at("sampleConfig");
    SampleControl samples(sampleConfigPath);
    for (const auto &varName : variables)
    {
        prepareHist(jsonConfig, varName, samples);
    }
    return 0;
}