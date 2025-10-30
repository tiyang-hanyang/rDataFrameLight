#include "SkimControl.h"
#include "Utility.h"

#include <iostream>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <memory>
#include <regex>
#include <csignal>
#include <unistd.h>
#include <atomic>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "ROOT/RSnapshotOptions.hxx"

////////////////////////////////////////////////// Setup configs

SkimControl::SkimControl(nlohmann::json configFile)
{
    this->readConfig(configFile);
}

SkimControl::SkimControl(const std::string &configPath)
{
    this->readConfig(configPath);
}

void SkimControl::readConfig(nlohmann::json origConfigFile)
{
    rdfWS_utility::JsonObject configFile(origConfigFile, "Skim JO config");
    // creating fileList would need era
    this->_ifMerging = configFile.at("merge").get<int>();
    this->_skimName = configFile.at("name").get<std::string>();
    this->_run = configFile.at("run").get<std::string>();
    this->_year = configFile.at("year").get<std::string>();
    this->_era = configFile.at("era").get<std::string>();

    // check if preliminary
    // by default is preliminary selection
    int preliminary(1);
    if (configFile.contains("preliminary"))
    {
        if (configFile.at("preliminary").get<int>() == 0) preliminary=0;
    }
    this->_isPreliminary = preliminary;

    // mkdir for output
    this->_outDir = configFile.at("outDir").get<std::string>();
    rdfWS_utility::creatingFolder("SkimControl", this->_outDir);

    // channel info
    this->_channels = configFile.at("datasets").get<std::vector<std::string>>();
    // by default, turn on all the channel skim
    for (const auto &channel : this->_channels)
    {
        this->_isOn.emplace(channel, 1);
    }
    // this->_isData = configFile.at("isData").get<std::map<std::string, int>>();
    this->_isData = configFile.at("isData").get<int>();
    this->_mcWeight = configFile.at("mcWeight").get<std::string>();

    // XS value should be in the json folder
    std::string XSConfigPath = configFile.at("XSConfig");
    if (XSConfigPath == "")
        this->_XSvalues = rdfWS_utility::readJson("SkimControl", "json/XS/" + this->_run + ".json");
    else
        this->_XSvalues = rdfWS_utility::readJson("SkimControl", XSConfigPath);
    // files
    std::string sampleConfigPath = configFile.at("sampleConfig");
    if (sampleConfigPath == "")
        this->_samples.emplace("json/samples/" + this->_era + ".json");
    else
        this->_samples.emplace(sampleConfigPath);
    // cuts
    std::vector<std::string> cutConfigList = configFile.at("cut");
    if (cutConfigList.size() > 0)
    {
        this->_skimCut = CutControl(cutConfigList[0]);
        for (int i = 1; i < cutConfigList.size(); i++)
        {
            this->_skimCut = this->_skimCut + CutControl(cutConfigList[i]);
        }
    }

    // branch list
    std::string branchPath = configFile.at("branchConfig");
    rdfWS_utility::JsonObject branchJson(rdfWS_utility::readJson("SkimControl", branchPath), "Branch Config");
    this->_branchList = branchJson.get<std::vector<std::string>>();
    // this->_branchList = configFile["branch"];
}

void SkimControl::readConfig(const std::string &configPath)
{
    auto jsonData = rdfWS_utility::readJson("SkimControl", configPath);
    this->readConfig(jsonData);
}

////////////////////////////////////////////////// Setup Data GoldenJson Lambda

void SkimControl::_createGoldenJsonFunc()
// ROOT::RDF::RNode SkimControl::applyGoldenJson(ROOT::RDF::RNode origData)
{
    rdfWS_utility::messageINFO("SkimControl", "Creating the golden json lambda function.");
    // check the golden json re-assignment
    if (this->_goldenJsonLambda)
    {
        rdfWS_utility::messageWARN("SkimControl", "The _goldenJsonLambda is already initialized. It is re-initialized again. Please make sure this is really what desired.");
    }

    // parsing the goldenJson
    std::string goldenJsonPath = "json/goldenJson/" + this->_year + "/goldenJson.json";
    // rdfWS_utility::messageINFO("SkimControl", "Loading golden json from: " + goldenJsonPath);
    auto goldenJson = rdfWS_utility::readJson("SkimControl", goldenJsonPath);
    std::map<std::string, std::vector<std::pair<int, int>>> goldenJsonList = goldenJson;

    // check the run number and then check the lumiblock
    this->_goldenJsonLambda =
        [goldenJsonList](ROOT::RDF::RNode origData)
    { return origData.Filter(
          [goldenJsonList](unsigned int run, unsigned int luminosityBlock)
          {
              if (goldenJsonList.find(std::to_string(run)) == goldenJsonList.end())
                  return false;
              auto goldenLumiBlk = goldenJsonList.at(std::to_string(run));
              int isGood(0);
              for (const auto &blkPair : goldenLumiBlk)
              {
                  if (luminosityBlock < blkPair.second && luminosityBlock > blkPair.first)
                  {
                      isGood = 1;
                      break;
                  }
              }
              if (!isGood)
                  return false;
              return true;
          },
          {"run", "luminosityBlock"}); };
}

void SkimControl::turnOn(const std::string &channels)
{
    bool found = false;

    std::regex pattern(channels);
    try
    {
        std::regex pattern(channels);

        for (auto &[key, value] : this->_isOn)
        {
            if (std::regex_match(key, pattern))
            {
                value = 1;
                found = true;
                rdfWS_utility::messageINFO("SkimControl", "Turned on channel: " + key);
            }
        }
    }
    catch (const std::regex_error &)
    {
        rdfWS_utility::messageWARN("SkimControl", "Invalid channel pattern: " + channels);
        return;
    }

    if (!found)
    {
        rdfWS_utility::messageWARN("SkimControl", "No channels matched: " + channels);
    }
}

void SkimControl::turnOff(const std::string &channels)
{
    bool found = false;

    std::regex pattern(channels);
    try
    {
        std::regex pattern(channels);

        for (auto &[key, value] : this->_isOn)
        {
            if (std::regex_match(key, pattern))
            {
                value = 0;
                found = true;
                rdfWS_utility::messageINFO("SkimControl", "Turned off channel: " + key);
            }
        }
    }
    catch (const std::regex_error &)
    {
        rdfWS_utility::messageWARN("SkimControl", "Invalid channel pattern: " + channels);
        return;
    }

    if (!found)
    {
        rdfWS_utility::messageWARN("SkimControl", "No channels matched: " + channels);
    }
}

////////////////////////////////////////////////// Exiting check gracefully

SkimControl *SkimControl::instance = nullptr;

void SkimControl::signalHandler(int signum)
{
    if (instance)
    {
        instance->stop_requested = true;
    }
}

////////////////////////////////////////////////// Skimming methods
// function for getting the total weight if exist
double SkimControl::_getTotalGenWeight(std::vector<std::string> fileLists)
{
    double totalGenWeight(0.0);
    for (const auto& fPath: fileLists)
    {
        std::unique_ptr<TFile> fTemp{ TFile::Open(fPath.c_str(), "READ") };
        if (!fTemp || fTemp->IsZombie()) 
        {
            rdfWS_utility::messageWARN("SkimControl", fPath+" not exist!");
            continue;
        }

        auto sumWeightHist = dynamic_cast<TH1D*>(fTemp->Get("genWeightSum"));
        if (!sumWeightHist)
        {
            rdfWS_utility::messageWARN("SkimControl", fPath+" does not have genWeightSum! Try manually add.");
            auto tTemp = dynamic_cast<TTree*>(fTemp->Get("Events"));

            ROOT::RDataFrame dfTemp(*tTemp);
            try {
                totalGenWeight += dfTemp.Sum("genWeight").GetValue();
            } catch (const std::exception& e) {
                rdfWS_utility::messageWARN("SkimControl", fPath+" cannot add genWeight, skip...");
                continue;
            }
        }
        else
            totalGenWeight += sumWeightHist->GetBinContent(1);
    }
    rdfWS_utility::messageINFO("SkimControl", "Total genWeight: " + std::to_string(totalGenWeight) );
    return totalGenWeight;
}

ROOT::RDF::RNode SkimControl::_preliminaryDeco(ROOT::RDF::RNode rndDS, int isData, const std::string& channel, double totalGenWeight) 
{
    // for preliminary skim, must applying golden json for data, and rescaling MC weight
    if (isData)
    {
        // apply golden json
        if (!(this->_goldenJsonLambda))
        {
            this->_createGoldenJsonFunc();
        }
        rndDS = this->_goldenJsonLambda.value()(rndDS);
    }
    else
    {
        // weight_XS = genWeight * (1 (lumi/fb-1) * 1000 * XS (pb)) / totalGenWeight
        double weightXSScale = 1000 * this->_XSvalues.at(channel) / totalGenWeight;
        std::ostringstream oss;
        oss << std::setprecision(15) << weightXSScale;
        std::string weightScaleStr = oss.str();
        rndDS = rndDS.Define("weight_XS", std::string("double(genWeight * ") + weightScaleStr + ")");
    }

    return rndDS;
}

std::vector<std::string> SkimControl::_getBranchArray(ROOT::RDF::RNode rndDS, int isData, int isPreliminary)
{
    // keep the branches in the config only and dump into files
    std::vector<std::string> branchArray;
    std::vector<std::string> originalBRs = rndDS.GetColumnNames();
    for (const auto &brName : this->_branchList)
    {
        if (std::find(originalBRs.begin(), originalBRs.end(), brName) == originalBRs.end())
            continue;
        branchArray.push_back(brName);
    }

    // MC need weight_XS
    if ((!isData) && isPreliminary) branchArray.push_back("weight_XS");

    return branchArray;
}

////////////////////////////////////////////////// Process the Skim running

void SkimControl::run()
{
    // enhance ctrl+C trigger sensitivity
    instance = this;
    signal(SIGINT, SkimControl::signalHandler);

    for (const auto &channel : this->_channels)
    {
        // detect stopping message, handling after this loop
        if (stop_requested)
        {
            rdfWS_utility::messageINFO("SkimControl", "Stop requested. Exiting before processing dataset " + channel);
            break;
        }

        // skip off channels 
        if (!this->_isOn[channel])
            continue;
        
        rdfWS_utility::messageINFO("SkimControl", "Starting processing channel " + channel);
        auto filePaths = this->_samples.value().getFiles(channel);
        if (filePaths.size() == 0)
            continue;
        
        // skim property
        auto isData = this->_isData;
        auto isPreSkim = this->_isPreliminary;

        // in case of MC, for preliminary skimming, need get the MC total genWeight
        double totalGenWeight(0.0);
        if ((!isData) && isPreSkim)
            totalGenWeight = _getTotalGenWeight(filePaths);

        // determine the outDir
        auto outDir = this->_outDir+"/"+ (isData?"data":"mc") +"/"+this->_era+"/"+channel+"/";
        rdfWS_utility::creatingFolder("SkimControl", outDir);

        // enable choose to do merging or not
        if (this->_ifMerging) 
        {
            // loading dataframe
            TChain *chDS = new TChain("Events");
            for (const auto &filePath : filePaths)
            {
                chDS->Add(filePath.c_str());
            }
            ROOT::RDataFrame rdfDS(*chDS);
            ROOT::RDF::RNode rndDS(rdfDS);

            // for preliminary skim, must applying golden json for data, and rescaling MC weight
            if (this->_isPreliminary)
                rndDS = this->_preliminaryDeco(rndDS, isData, channel, totalGenWeight);

            // apply the filter
            rndDS = this->_skimCut.applyCut(rndDS);
            if (rndDS.Count().GetValue() == 0)
                continue;

            // keep the branches in the config only and dump into files
            auto branchArray = this->_getBranchArray(rndDS, isData, isPreSkim);

            // output
            ROOT::RDF::RSnapshotOptions SSoption;
            SSoption.fCompressionLevel = 9;
            std::string outputPath = outDir + "/" + channel + "_" + this->_skimName + ".root";
            rndDS.Snapshot("Events", outputPath, branchArray, SSoption);

            // collect memory
            delete chDS;
        }
        else
        {
            for (const auto &filePath: filePaths)
            {
                // check if alreayd processed first
                // file naming structure
                // original: data or mc /era/channel/ (NANOAOD/MINIAOD) or condition/ runNumber / sampleName
                // after Rochester corr, uniformed: data or mc /era/channel/ (runNumber+sample)
                auto outSampleName = filePath.substr(filePath.rfind("/")+1);
                if (!outSampleName.find("Rcorr")) 
                {
                    auto dirPart = filePath.substr(0, filePath.rfind("/"));
                    auto runNumber =  dirPart.substr(dirPart.rfind("/")+1);
                    outSampleName = runNumber+"-"+outSampleName;
                }
                std::string outputPath = outDir + "/" + outSampleName;

                if (std::filesystem::exists(outputPath)) 
                {
                    rdfWS_utility::messageINFO("SkimControl", outputPath + " already exist, skip");
                    continue;
                }

                // loading dataframe
                ROOT::RDataFrame rdfDS("Events", filePath.c_str());
                ROOT::RDF::RNode rndDS(rdfDS);

                // for preliminary skim, must applying golden json for data, and rescaling MC weight
                if (this->_isPreliminary)
                    rndDS = this->_preliminaryDeco(rndDS, isData, channel, totalGenWeight);

                // apply the filter
                rndDS = this->_skimCut.applyCut(rndDS);
                if (rndDS.Count().GetValue() == 0)
                    continue;

                // keep the branches in the config only and dump into files
                auto branchArray = this->_getBranchArray(rndDS, isData, isPreSkim);

                // output
                ROOT::RDF::RSnapshotOptions SSoption;
                SSoption.fCompressionLevel = 9;
                rndDS.Snapshot("Events", outputPath, branchArray, SSoption);
            }
        }
    }
}