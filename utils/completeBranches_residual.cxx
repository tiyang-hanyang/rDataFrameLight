// make_muon_friend.C
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TFile.h>
#include <TTree.h>
#include <string>
#include <map>
#include <filesystem>
#include <algorithm>
#include <fstream>

#include "SampleControl.h"
#include "Utility.h"

std::map<std::string, std::string> getTarget(std::vector<std::string> skimmedFilePaths, const std::string& channel)
{
    std::map<std::string, std::string> outMap;
    std::filesystem::create_directory("/home/tiyang/public/br_complete_residual/");
    for (const auto &path : skimmedFilePaths)
    {
        if (path.find(channel) == std::string::npos)
        {
            rdfWS_utility::messageWARN("completeBranches", "The file path "+path+" does not contain the channel "+channel+", skip...");
            continue;
        }
        std::string folder = path.substr(path.find(channel));
        folder = folder.substr(0, folder.rfind("/"));
        std::string toGen = "/home/tiyang/public/br_complete_residual/"+folder;
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
                       const std::map<std::string, std::string> &targetFileMap,
                       std::size_t &written,
                       std::size_t batchLimit)
{
    for (const auto &kv : targetFileMap)
    {
        const auto &skimmedFile = kv.first;
        const auto &outFile = kv.second;
        ofs << skimmedFile << " " << outFile << "\n";
        ++written;
        if (batchLimit > 0 && written >= batchLimit)
            return true;
    }
    return true;
}

int main(int argc, char *argv[])
{
    // samples
    std::string skimmedJson = "/home/tiyang/public/rDataFrameLight_update/source/json/samples/Dimuon_NanoAOD/RunIII2024Summer24NanoAODv15_forCorr_temp.json";

    std::string jobOutBase = "/home/tiyang/public/br_complete_residual/jobs";
    std::size_t batchLimit = 1500;

    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        if (arg == "--skimmedJson" && i + 1 < argc)
        {
            skimmedJson = argv[++i];
        }
        else if (arg == "--jobOut" && i + 1 < argc)
        {
            jobOutBase = argv[++i];
        }
        else if (arg == "--batchLimit" && i + 1 < argc)
        {
            batchLimit = static_cast<std::size_t>(std::stoul(argv[++i]));
        }
        else
        {
            rdfWS_utility::messageWARN("completeBranches",
                                       "Unknown or incomplete argument: " + arg);
            return 1;
        }
    }

    SampleControl skimmedSamplesArg(skimmedJson);

    // generate folder for jobOutBase.
    const std::filesystem::path outPath(jobOutBase);
    if (!outPath.parent_path().empty() && !std::filesystem::exists(outPath.parent_path()))
    {
        std::filesystem::create_directories(outPath.parent_path());
    }

    auto allChannels = skimmedSamplesArg.getAllChannels();
    std::size_t batchIndex = 0;
    std::size_t writtenInBatch = 0;
    std::ofstream ofs;
    auto openBatch = [&]() -> bool {
        if (ofs.is_open())
            ofs.close();
        const std::string jobOut = jobOutBase + "_" + std::to_string(batchIndex) + ".txt";
        ofs.open(jobOut);
        if (!ofs.is_open())
        {
            rdfWS_utility::messageERROR("completeBranches", "Failed to open job file: " + jobOut);
            return false;
        }
        return true;
    };
    if (!openBatch())
        return 1;

    for (const auto &channel : allChannels)
    {
        auto skimmedFilePaths = skimmedSamplesArg.getFiles(channel);
        auto targetFileMap = getTarget(skimmedFilePaths, channel);
        while (true)
        {
            if (!appendJobs(ofs, targetFileMap, writtenInBatch, batchLimit))
                return 1;
            if (batchLimit == 0 || writtenInBatch < batchLimit)
                break;
            ++batchIndex;
            writtenInBatch = 0;
            if (!openBatch())
                return 1;
        }
    }
    return 0;
}
