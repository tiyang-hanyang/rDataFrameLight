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

////////////////////////////////////////////////// Process the Skim running

void SkimControl::run()
{
    // add a ctrl+C detector
    instance = this;
    signal(SIGINT, SkimControl::signalHandler);

    // need to use chain, as skim requires reweighting based
    for (const auto &channel : this->_channels)
    {
        // detect stopping message, handling after this loop
        if (stop_requested)
        {
            rdfWS_utility::messageINFO("SkimControl", "Stop requested. Exiting before processing dataset " + channel);
            break;
        }

        if (!this->_isOn[channel])
            continue;
        rdfWS_utility::messageINFO("SkimControl", "Starting processing channel " + channel);
        auto filePaths = this->_samples.value().getFiles(channel);
        if (filePaths.size() == 0)
            continue;
        // if is data
        auto isData = this->_isData;

        // loading dataframe
        TChain *chDS = new TChain("Events");
        for (const auto &filePath : filePaths)
        {
            chDS->Add(filePath.c_str());
        }
        ROOT::RDataFrame rdfDS(*chDS);
        ROOT::RDF::RNode rndDS(rdfDS);

        // determine output Dir for this channel
        // non preliminary not need to have the folder spliting, as not doing MC reweighting or data golden json match
        auto outDir = this->_outDir;
        if (this->_isPreliminary)
        {
            if (isData)
            {
                std::string dirPart = filePaths[0];
                if (dirPart.find("/data/") == std::string::npos)
                    rdfWS_utility::messageERROR("SkimControl", "The data path does not have '/data/' in the path, please have a check");
                dirPart = dirPart.substr(dirPart.find("/data/"));
                dirPart = dirPart.substr(0, dirPart.rfind("/"));
                // 2024 storage structure is different, need special treatment
                // a typical 2023 sample: /data1/common/NanoAOD/data/Run2023C/Muon0/NANOAOD/22Sep2023_v1-v1/2530000/030eead5-93f9-405c-863c-a62244712e91.root
                // 2024 case: /data1/common/NanoAOD/data/Run2024C/Muon0/NANOAOD/PromptReco-v1/000/379/415/00000/80b3773e-c489-435b-a25b-027ff8e64f64.root
                if (this->_year == "2024")
                {
                    dirPart = dirPart.substr(0, dirPart.rfind("/"));
                    dirPart = dirPart.substr(0, dirPart.rfind("/"));
                    dirPart = dirPart.substr(0, dirPart.rfind("/"));
                }
                dirPart = dirPart.substr(0, dirPart.rfind("/") + 1);
                outDir += dirPart;
            }
            else
            {
                std::string dirPart = filePaths[0];
                if (dirPart.find("/mc/") == std::string::npos)
                    rdfWS_utility::messageERROR("SkimControl", "The data path does not have '/mc/' in the path, please have a check");
                dirPart = dirPart.substr(dirPart.find("/mc/"));
                dirPart = dirPart.substr(0, dirPart.rfind("/"));
                dirPart = dirPart.substr(0, dirPart.rfind("/") + 1);
                outDir += dirPart;
            }

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
                // need to compute weight of XS before taking the filter for MC
                double totalGenWeight = rndDS.Sum("genWeight").GetValue();
std::cout << "total sum gen weight: " << totalGenWeight << std::endl;
                // store weight to XS per lumi to enable scaling
                // weight_XS = genWeight * (1 (lumi/fb-1) * 1000 * XS (pb)) / totalGenWeight
                double weightXSScale = 1000 * this->_XSvalues.at(channel) / totalGenWeight;

//std::cout << "weight XS scale number: " << std::to_string(weightXSScale) << std::endl;
//                rndDS = rndDS.Define("weight_XS", std::string("(double(genWeight) * double(") + std::to_string(weightXSScale) + "))");
                std::ostringstream oss;
                oss << std::setprecision(15) << weightXSScale;
                std::string weightScaleStr = oss.str();
std::cout << "weight XS scale number: " << weightScaleStr << std::endl;
                rndDS = rndDS.Define("weight_XS", std::string("(genWeight * ") + weightScaleStr + ")");

std::cout << "Now, tell me what is the weight_XS after the scaling: " << rndDS.Sum("weight_XS").GetValue() << std::endl;
            }
        }

        rdfWS_utility::creatingFolder("SkimControl", outDir + "/" + this->_era + "/");
        std::string outputPath = outDir + "/" + this->_era + "/" + channel + "_" + this->_skimName + ".root";

        // apply the filter
        rndDS = this->_skimCut.applyCut(rndDS);
        if (rndDS.Count().GetValue() == 0)
            continue;

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
        // if (!(this->_isData.at(channel)))
        if (!isData)
            branchArray.push_back("weight_XS");

        // output
        ROOT::RDF::RSnapshotOptions SSoption;
        SSoption.fCompressionLevel = 9;
        rndDS.Snapshot("Events", outputPath, branchArray, SSoption);
    }
}
