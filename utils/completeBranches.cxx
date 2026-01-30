// make_muon_friend.C
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TFile.h>
#include <TTree.h>
#include <unordered_set>
#include <string>
#include <map>
#include <filesystem>
#include <algorithm>
#include <fstream>

#include "SampleControl.h"
#include "Utility.h"

// function to accumulate the in-out matching map
std::map<std::string, std::string> getSource(const std::vector<std::string> &skimmedFilePaths, const std::vector<std::string> &rawFilePaths)
{
    std::map<std::string, std::string> rawMap;
    for (const auto &path : skimmedFilePaths)
    {
        // e.g. /data2/common/skimmed_NanoAOD/skim_0126_medium_dimuon/mc/RunIII2024Summer24NanoAODv15/DY2Mu/100000-a1c5529b-7c90-4595-8f50-c74d9e596494.root
        std::string fileName = path.substr(path.rfind("/") + 1);
        std::string origFileName = fileName.substr(fileName.find("-") + 1);
        for (const auto &rawPath : rawFilePaths)
        {
            // e.g. "/data2/common/NanoAOD/mc/v15/RunIII2024Summer24NanoAODv15/DYto2Mu-2Jets_Bin-MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v6/100000/a1c5529b-7c90-4595-8f50-c74d9e596494.root"
            if (rawPath.find(origFileName) != std::string::npos)
            {
                rawMap.emplace(path, rawPath);
                continue;
            }
        }
    }
    return rawMap;
}

std::map<std::string, std::string> getTarget(std::vector<std::string> skimmedFilePaths, const std::string& channel)
{
    std::map<std::string, std::string> outMap;
    std::filesystem::create_directory("/home/tiyang/public/br_complete/");
    for (const auto &path : skimmedFilePaths)
    {
        if (path.find(channel) == std::string::npos)
        {
            rdfWS_utility::messageWARN("completeBranches", "The file path "+path+" does not contain the channel "+channel+", skip...");
            continue;
        }
        std::string folder = path.substr(path.find(channel));
        folder = folder.substr(0, folder.rfind("/"));
        std::string toGen = "/home/tiyang/public/br_complete/"+folder;
        if (!std::filesystem::exists(toGen))
        {
            std::filesystem::create_directories(toGen);
        }
        std::string fileName = path.substr(path.rfind("/") + 1);
        std::string outPath = toGen + "/" + fileName;
        outMap.emplace(path, outPath);
    }
    return outMap;
}

static bool appendJobs(std::ofstream &ofs,
                       const std::map<std::string, std::string> &rawFileMap,
                       const std::map<std::string, std::string> &targetFileMap,
                       const std::string &branchJson)
{
    for (const auto &kv : rawFileMap)
    {
        const auto &skimmedFile = kv.first;
        const auto &rawFile = kv.second;
        auto it = targetFileMap.find(skimmedFile);
        if (it == targetFileMap.end())
        {
            rdfWS_utility::messageWARN("completeBranches", "No output path for skim file: " + skimmedFile);
            continue;
        }
        const auto &outFile = it->second;
        ofs << skimmedFile << " " << rawFile << " " << outFile << " " << branchJson << "\n";
    }
    return true;
}

int main(int argc, char *argv[])
{
    // samples
    std::string skimmedJson = "/home/tiyang/public/rDataFrameLight_update/source/json/samples/Dimuon_NanoAOD/RunIII2024Summer24NanoAODv15_forCorr_temp.json";
    std::string rawJson = "/home/tiyang/public/rDataFrameLight_update/source/json/samples/rawNanoAOD/RunIII2024Summer24NanoAODv15.json";

    // branches
    std::string branchJson = "/home/tiyang/public/rDataFrameLight_update/source/json/branches/MuonMetricBranchNew.json";
    std::string jobOut = "/home/tiyang/public/br_complete/jobs.txt";

    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        if (arg == "--skimmedJson" && i + 1 < argc)
        {
            skimmedJson = argv[++i];
        }
        else if (arg == "--rawJson" && i + 1 < argc)
        {
            rawJson = argv[++i];
        }
        else if (arg == "--branches" && i + 1 < argc)
        {
            branchJson = argv[++i];
        }
        else if (arg == "--jobOut" && i + 1 < argc)
        {
            jobOut = argv[++i];
        }
        else
        {
            rdfWS_utility::messageWARN("completeBranches",
                                       "Unknown or incomplete argument: " + arg);
            return 1;
        }
    }

    SampleControl skimmedSamplesArg(skimmedJson);
    SampleControl rawSamplesArg(rawJson);

    // generate folder for jobOut.
    const std::filesystem::path outPath(jobOut);
    if (!outPath.parent_path().empty() && !std::filesystem::exists(outPath.parent_path()))
    {
        std::filesystem::create_directories(outPath.parent_path());
    }

    std::ofstream ofs(jobOut);
    if (!ofs.is_open())
    {
        rdfWS_utility::messageERROR("completeBranches", "Failed to open job file: " + jobOut);
        return 1;
    }

    auto allChannels = skimmedSamplesArg.getAllChannels();
    for (const auto &channel : allChannels)
    {
        auto skimmedFilePaths = skimmedSamplesArg.getFiles(channel);
        auto rawFilePaths = rawSamplesArg.getFiles(channel);
        auto rawFileMap = getSource(skimmedFilePaths, rawFilePaths);
        auto targetFileMap = getTarget(skimmedFilePaths, channel);
        if (!appendJobs(ofs, rawFileMap, targetFileMap, branchJson))
            return 1;
    }
    return 0;
}
