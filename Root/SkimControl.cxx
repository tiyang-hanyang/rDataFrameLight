#include "SkimControl.h"
#include "Utility.h"

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <memory>
#include <regex>
#include <csignal>
#include <unistd.h>
#include <atomic>
#include <unordered_map>
#include <algorithm>
#include <cctype>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "ROOT/RSnapshotOptions.hxx"
#include "ROOT/RDFHelpers.hxx"

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

    // job type check, better to have to avoid confusion
    if (configFile.contains("jobType"))
    {
        std::string jobType = configFile.at("jobType").get<std::string>();
        if (jobType != "skim") 
        {
            rdfWS_utility::messageERROR("SkimControl.cxx", "The jobType of your config is not skimming! Please check again. Running ceases...");
            exit(1);
        }
    }

    this->_ifMerging = configFile.at("merge").get<int>();
    this->_year = configFile.at("year").get<std::string>();
    this->_era = configFile.at("era").get<std::string>();
    this->_skimName = this->_era + "_skimmmed";
    if (configFile.contains("name"))
    {
       this->_skimName = configFile.at("name").get<std::string>();
    }

    // by default do preliminary selection
    // mc sum genWeight from Events tree
    // data need test the golden json
    int preliminary(1);
    if (configFile.contains("preliminary"))
    {
        if (configFile.at("preliminary").get<int>() == 0) preliminary=0;
    }
    this->_isPreliminary = preliminary;

    // set maximum files for skimming
    int maxFiles(0);
    if (configFile.contains("maxFiles"))
    {
        maxFiles = configFile.at("maxFiles").get<int>();
        if (maxFiles < 0) maxFiles = 0;
    }
    this->_maxFilesPerChannel = maxFiles;

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

    // files
    std::string sampleConfigPath = configFile.at("sampleConfig");
    if (sampleConfigPath == "")
        this->_samples.emplace("json/samples/" + this->_era + ".json");
    else
        this->_samples.emplace(sampleConfigPath);

    // prepare all the cuts
    std::vector<std::string> cutConfigList = configFile.at("cutConfig");
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
}

void SkimControl::readConfig(const std::string &configPath)
{
    auto jsonData = rdfWS_utility::readJson("SkimControl", configPath);
    this->readConfig(jsonData);
}

////////////////////////////////////////////////// Setup Data GoldenJson Lambda
void SkimControl::_createGoldenJsonFunc()
{
    rdfWS_utility::messageINFO("SkimControl", "Creating the golden json lambda function.");
    // check the golden json re-assignment
    if (this->_goldenJsonLambda)
    {
        rdfWS_utility::messageWARN("SkimControl", "The _goldenJsonLambda is already initialized. It is re-initialized again. Please make sure this is really what desired.");
    }

    // parsing the goldenJson, cast from string to int
    std::string goldenJsonPath = "json/goldenJson/" + this->_year + "/goldenJson.json";
    // rdfWS_utility::messageINFO("SkimControl", "Loading golden json from: " + goldenJsonPath);
    auto goldenJson = rdfWS_utility::readJson("SkimControl", goldenJsonPath);
    std::map<std::string, std::vector<std::pair<unsigned int, unsigned int>>> goldenJsonRaw = goldenJson;
    std::unordered_map<unsigned int, std::vector<std::pair<unsigned int,unsigned int>>> goldenJsonList;
    goldenJsonList.reserve(goldenJsonRaw.size());

    auto is_digits = [](const std::string& s) {
        return !s.empty() && std::all_of(s.begin(), s.end(), [](unsigned char c) { return std::isdigit(c) != 0; });
    };

    for (const auto& [k,v] : goldenJsonRaw) 
    {
        if (!is_digits(k)) {
            rdfWS_utility::messageWARN("SkimControl", "Invalid run key: " + k + ", skip.");
            continue;
        }
        unsigned int run = std::stoul(k);
        goldenJsonList.emplace(run, v);
    }

    // golden json function definition
    this->_goldenJsonLambda =
        [goldenJsonList](ROOT::RDF::RNode origData)
    { return origData.Filter(
        [goldenJsonList](unsigned int run, unsigned int luminosityBlock)
        {
            if (goldenJsonList.find(run) == goldenJsonList.end())
                return false;
            auto goldenLumiBlk = goldenJsonList.at(run);
            int isGood(0);
            for (const auto &blkPair : goldenLumiBlk)
            {
                if (luminosityBlock <= blkPair.second && luminosityBlock >= blkPair.first)
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

////////////////////////////////////////////////// Enable exiting gracefully

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
            if (!fTemp->GetListOfKeys()->FindObject("Events")) continue;
            auto tTemp = dynamic_cast<TTree*>(fTemp->Get("Events"));
            if (!tTemp) 
            {
                rdfWS_utility::messageWARN("SkimControl", fPath+" Events is not a TTree, skip...");
                continue;
            }

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

ROOT::RDF::RNode SkimControl::_preliminaryDeco(ROOT::RDF::RNode rndDS, const std::string& channel, double totalGenWeight) 
{
    // apply golden json
    if (!(this->_goldenJsonLambda))
    {
        this->_createGoldenJsonFunc();
    }
    rndDS = this->_goldenJsonLambda.value()(rndDS);

    return rndDS;
}

std::vector<std::string> SkimControl::_getBranchArray(ROOT::RDF::RNode rndDS, int isPreliminary)
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

    return branchArray;
}

////////////////////////////////////////////////// Process the Skim running

void SkimControl::run()
{
    // enhance ctrl+C trigger sensitivity
    instance = this;
    signal(SIGINT, SkimControl::signalHandler);
    signal(SIGTERM, SkimControl::signalHandler);


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
        // if ((!isData) && isPreSkim)
        //     totalGenWeight = _getTotalGenWeight(filePaths);

        // determine the outDir
        auto outDir = this->_outDir+"/"+ (isData?"data":"mc") +"/"+this->_era+"/"+channel+"/";
        rdfWS_utility::creatingFolder("SkimControl", outDir);

        // enable choose to do merging or not
        if (this->_ifMerging) 
        {
            // collect total genWeight for merged output (MC only, preliminary skim)
            double totalGenWeightMerged = 0.0;
            if (!isData && isPreSkim)
                totalGenWeightMerged = _getTotalGenWeight(filePaths);

            // loading dataframe with chain safety check
            TChain *chDS = new TChain("Events");
            for (const auto &filePath : filePaths)
            {
                std::unique_ptr<TFile> fTemp{ TFile::Open(filePath.c_str(), "READ") };
                if (!fTemp || fTemp->IsZombie())
                {
                    rdfWS_utility::messageWARN("SkimControl", filePath + " not exist!");
                    continue;
                }
                // for the case of non preliminary skimming, I must having the genWeightSum
                if (!isData && !isPreSkim)
                {
                    if (!fTemp->GetListOfKeys()->FindObject("genWeightSum"))
                    {
                        rdfWS_utility::messageERROR("SkimControl", filePath + " has no genWeightSum for non-preliminary skimming. This is very strange, stop running, please check!");
                        return;
                    }
                    auto tempsumhist = dynamic_cast<TH1D*>(fTemp->Get("genWeightSum"));
                    if (!tempsumhist)
                    {
                        rdfWS_utility::messageERROR("SkimControl", filePath + " has not valid genWeightSum for non-preliminary skimming. This is very strange, stop running, please check!");
                        return;
                    }
                    totalGenWeightMerged += tempsumhist->Integral();
                }

                if (!fTemp->GetListOfKeys()->FindObject("Events"))
                {
                    rdfWS_utility::messageWARN("SkimControl", filePath + " has no Events tree, skip chaining.");
                    continue;
                }
                chDS->Add(filePath.c_str());
            }
            ROOT::RDataFrame rdfDS(*chDS);
            ROOT::RDF::RNode rndDS(rdfDS);

            // for preliminary skim, must applying golden json for data
            if (this->_isPreliminary && isData)
                rndDS = this->_preliminaryDeco(rndDS, channel, totalGenWeightMerged);

            // apply the filter
            rndDS = this->_skimCut.applyCut(rndDS);

            // keep the branches in the config only and dump into files
            auto branchArray = this->_getBranchArray(rndDS, isPreSkim);

            // output
            ROOT::RDF::RSnapshotOptions SSoption;
            SSoption.fCompressionLevel = 6;
            std::string outputPath = outDir + "/" + channel + "_" + this->_skimName + ".root";
            rndDS.Snapshot("Events", outputPath, branchArray, SSoption);

            if (!isData)
            {
                TFile* fRC = new TFile(outputPath.c_str(), "UPDATE");
                auto h_sumw = new TH1D("genWeightSum", "sum of genWeight (this file)", 1, 0.0, 1.0);
                h_sumw->SetBinContent(1, totalGenWeightMerged);
                h_sumw->GetYaxis()->SetTitle("sum(genWeight)");
                h_sumw->Write("");
                fRC->Close();
                delete fRC;
            }

            // collect memory
            delete chDS;
        }
        else
        {
            std::vector<std::pair<std::string, std::string>> jobs;

            for (const auto &filePath: filePaths)
            {
                // file naming structure
                // original: data or mc /era/channel/ (NANOAOD/MINIAOD) or condition/ runNumber / sampleName
                // after Rochester corr, uniformed: data or mc /era/channel/ (runNumber+sample)
                auto outSampleName = filePath.substr(filePath.rfind("/")+1);
                auto dirPart = filePath.substr(0, filePath.rfind("/"));
                auto runNumber =  dirPart.substr(dirPart.rfind("/")+1);
                outSampleName = runNumber+"-"+outSampleName;

                std::string outputPath = outDir + "/" + outSampleName;

                if (std::filesystem::exists(outputPath)) 
                {
                    rdfWS_utility::messageINFO("SkimControl", outputPath + " already exist, skip");
                    continue;
                }
                jobs.emplace_back(filePath, outputPath);
            }

            if (this->_maxFilesPerChannel>0 && jobs.size()>this->_maxFilesPerChannel) 
            {
                jobs.erase(jobs.begin()+this->_maxFilesPerChannel, jobs.end());
            }

            auto totalJobs = jobs.size();
            if (totalJobs==0) continue;
            int batchSize=8; // hard coding now

            // batchSize=1;
            // doing parallel skimming, after each batch, need to write back the weight back
            if (batchSize>1)
            {
                int batches = (totalJobs-1)/batchSize+1;

                // running inside a batch
                for (int i=0; i < batches; i++)
                {
                    rdfWS_utility::messageINFO("SkimControl", "Executing batch "+std::to_string(i));

                    if (stop_requested)
                    {
                        rdfWS_utility::messageINFO("SkimControl", "Stop requested. Exiting before processing dataset " + channel);
                        break;
                    }

                    std::vector<std::unique_ptr<ROOT::RDataFrame>> dfs;

                    using SnapRet_t = ROOT::RDF::RResultPtr<ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager>>;
                    using SumRet_t  = ROOT::RDF::RResultPtr<float>;
                    std::vector<SnapRet_t> snapshots; 
                    std::vector<SumRet_t>  sumws;

                    std::vector<std::string> outPaths; 

                    int nInBatch = batchSize;
                    if (i==batches-1) nInBatch = totalJobs-batchSize*i;

                    dfs.reserve(nInBatch);
                    snapshots.reserve(nInBatch);
                    sumws.reserve(nInBatch);
                    outPaths.reserve(nInBatch);

                    ROOT::RDF::RSnapshotOptions SSoption;
                    SSoption.fCompressionLevel = 4;
                    SSoption.fLazy = true;

                    for (int jobId = batchSize*i; jobId < batchSize*i+nInBatch; jobId++)
                    {
                        auto filePath = jobs[jobId].first;
                        auto outputPath = jobs[jobId].second;

                        // loading dataframe
                        // ROOT::RDataFrame rdfDS("Events", filePath.c_str());
                        rdfWS_utility::messageINFO("SkimControl", "Adding "+filePath);
                        dfs.emplace_back(std::make_unique<ROOT::RDataFrame>("Events", filePath.c_str()));
                        auto& rdfDS = *dfs.back();
                        ROOT::RDF::RNode rndDS(rdfDS);

                        // note could only for data
                        if (!isData)
                            sumws.emplace_back(rndDS.Sum<float>("genWeight"));
                        else
                            sumws.emplace_back(rndDS.Define("dummyWeight", "1.0f").Sum<float>("dummyWeight"));
                            // sumws.emplace_back(rndDS.Count());

                        // for preliminary skim, must applying golden json for data, and rescaling MC weight
                        if (this->_isPreliminary)
                            rndDS = this->_preliminaryDeco(rndDS, channel, totalGenWeight);

                        // apply the filter
                        rndDS = this->_skimCut.applyCut(rndDS);

                        // keep the branches in the config only and dump into files
                        auto branchArray = this->_getBranchArray(rndDS, isPreSkim);

                        // output
                        auto snap = rndDS.Snapshot("Events", outputPath, branchArray, SSoption);
                        outPaths.emplace_back(outputPath);
                        snapshots.push_back(snap);
                    }
                    std::vector<ROOT::RDF::RResultHandle> handles;
                    handles.reserve(snapshots.size() + sumws.size());
                    for (auto &s : snapshots) handles.emplace_back(ROOT::RDF::RResultHandle(s));
                    for (auto &w : sumws)     handles.emplace_back(ROOT::RDF::RResultHandle(w));
                    unsigned ran = ROOT::RDF::RunGraphs(handles);

                    // after a batch running, padding the total weight
                    for (int idx_inB = 0; idx_inB < nInBatch; idx_inB++)
                    {
                        auto outputPath = outPaths[idx_inB];
                        const double totalWeight = sumws[idx_inB].GetValue();

                        TFile* fRC = new TFile(outputPath.c_str(), "UPDATE");
                        auto h_sumw = new TH1D("genWeightSum", "sum of genWeight (this file)", 1, 0.0, 1.0);
                        h_sumw->SetBinContent(1, totalWeight);
                        h_sumw->GetYaxis()->SetTitle("sum(genWeight)");
                        h_sumw->Write("");
                        fRC->Close();
                        delete fRC;
                    }
                }
            }
            else
            {
                rdfWS_utility::messageINFO("SkimControl", "Executing non batch ");

                std::vector<std::unique_ptr<ROOT::RDataFrame>> dfs;

                using SnapRet_t = ROOT::RDF::RResultPtr<ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager>>;
                std::vector<SnapRet_t> snapshots; 

                ROOT::RDF::RSnapshotOptions SSoption;
                SSoption.fCompressionLevel = 4;
                SSoption.fLazy = false;

                for (int jobId = 0; jobId < totalJobs; jobId++)
                {
                    auto filePath = jobs[jobId].first;
                    auto outputPath = jobs[jobId].second;

                    if (stop_requested)
                    {
                        rdfWS_utility::messageINFO("SkimControl", "Stop requested. Exiting before processing dataset " + channel);
                        break;
                    }

                    // loading dataframe
                    // ROOT::RDataFrame rdfDS("Events", filePath.c_str());
                    rdfWS_utility::messageINFO("SkimControl", "Adding "+filePath);
                    dfs.emplace_back(std::make_unique<ROOT::RDataFrame>("Events", filePath.c_str()));
                    auto& rdfDS = *dfs.back();
                    ROOT::RDF::RNode rndDS(rdfDS);

                    // for preliminary skim, must applying golden json for data, and rescaling MC weight
                    if (this->_isPreliminary)
                        rndDS = this->_preliminaryDeco(rndDS, channel, totalGenWeight);

                    // apply the filter
                    rndDS = this->_skimCut.applyCut(rndDS);

                    // keep the branches in the config only and dump into files
                    auto branchArray = this->_getBranchArray(rndDS, isPreSkim);

                    // output
                    rndDS.Snapshot("Events", outputPath, branchArray, SSoption);
                }
            }
        }
    }
}
