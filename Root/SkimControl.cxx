#include "SkimControl.h"
#include "Utility.h"

#include <iostream>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <memory>
#include <regex>

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

void SkimControl::readConfig(nlohmann::json configFile)
{
    // doing the skim need: skimName, channelList, isData, XS_values
    // creating fileList would need era
    this->_skimName = configFile["name"];
    this->_run = configFile["run"];
    this->_year = configFile["year"];
    this->_era = configFile["era"];

    // mkdir for output
    this->_outDir = configFile["outDir"];
    rdfWS_utility::creatingFolder("[SkimControl]", this->_outDir);

    // channel info
    this->_channels = configFile["datasets"];
    // by default, turn on all the channel skim
    for (const auto &channel : this->_channels)
    {
        this->_isOn.emplace(channel, 1);
    }
    this->_isData = configFile["isData"];
    this->_mcWeight = configFile["mcWeight"];

    // XS value should be in the json folder
    if (configFile["XSConfig"] == "")
        this->_XSvalues = rdfWS_utility::readJson("SkimControl", "json/XS/" + this->_run + ".json");
    else
        this->_XSvalues = rdfWS_utility::readJson("SkimControl", configFile["XSConfig"]);
    // files
    if (configFile["sampleConfig"] == "")
        this->_samples.emplace("json/samples/" + this->_era + ".json");
    else
        this->_samples.emplace(configFile["sampleConfig"]);
    // cuts
    std::vector<std::string> cutConfigList = configFile["cut"];
    if (cutConfigList.size() > 0)
    {
        this->_skimCut = CutControl(cutConfigList[0]);
        for (int i = 1; i < cutConfigList.size(); i++)
        {
            this->_skimCut = this->_skimCut + CutControl(cutConfigList[i]);
        }
    }

    // branch list
    std::string branchConfig = configFile["branchConfig"];
    this->_branchList = nlohmann::json(branchConfig);
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

////////////////////////////////////////////////// Process the Skim running

void SkimControl::run()
{
    // need to use chain, as skim requires reweighting based
    for (const auto &channel : this->_channels)
    {
        if (!this->_isOn[channel])
            continue;
        rdfWS_utility::messageINFO("SkimControl", "Starting processing channel " + channel);
        auto filePaths = this->_samples.value().getFiles(channel);
        if (filePaths.size() == 0)
            continue;
        // if is data
        auto isData = this->_isData.at(channel);

        // loading dataframe
        TChain *chDS = new TChain("Events");
        for (const auto &filePath : filePaths)
        {
            chDS->Add(filePath.c_str());
        }
        ROOT::RDataFrame rdfDS(*chDS);
        ROOT::RDF::RNode rndDS(rdfDS);

        // determine output Dir for this channel
        auto outDir = this->_outDir;
        if (isData)
        {
            std::string dirPart = filePaths[0];
            dirPart = dirPart.substr(dirPart.find("/data/"));
            dirPart = dirPart.substr(0, dirPart.rfind("/"));
            dirPart = dirPart.substr(0, dirPart.rfind("/") + 1);
            outDir += dirPart;
        }
        else
        {
            std::string dirPart = filePaths[0];
            dirPart = dirPart.substr(dirPart.find("/mc/"));
            dirPart = dirPart.substr(0, dirPart.rfind("/"));
            dirPart = dirPart.substr(0, dirPart.rfind("/") + 1);
            outDir += dirPart;
        }
        rdfWS_utility::creatingFolder("SkimControl", outDir);
        std::string outputPath = outDir + "/" + channel + "_" + this->_skimName + ".root";

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
            float totalGenWeight = rndDS.Sum("genWeight").GetValue();
            // store weight to XS per lumi to enable scaling
            // weight_XS = genWeight * (1 (lumi/fb-1) * 1000 * XS (pb)) / totalGenWeight
            float weightXSScale = 1000 * this->_XSvalues.at(channel) / totalGenWeight;
            rndDS = rndDS.Define("weight_XS", std::string("(genWeight * ") + std::to_string(weightXSScale) + ")");
        }

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
        if (!(this->_isData.at(channel)))
            branchArray.push_back("weight_XS");

        // output
        ROOT::RDF::RSnapshotOptions SSoption;
        SSoption.fCompressionLevel = 9;
        rndDS.Snapshot("Events", outputPath, branchArray, SSoption);
    }
}