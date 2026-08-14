#include "Utility.h"
#include "SampleControl.h"

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "TFile.h"

namespace
{
std::string buildOutputPath(const std::string &outDir,
                            int isData,
                            const std::string &era,
                            const std::string &channel,
                            const std::string &filePath)
{
    auto outBase = outDir + "/" + (isData ? "data" : "mc") + "/" + era + "/" + channel + "/";
    auto outSampleName = filePath.substr(filePath.rfind("/") + 1);
    auto dirPart = filePath.substr(0, filePath.rfind("/"));
    auto runNumber = dirPart.substr(dirPart.rfind("/") + 1);
    outSampleName = runNumber + "-" + outSampleName;
    return outBase + "/" + outSampleName;
}

static bool appendJobs(std::ofstream &ofs,
                       const std::vector<std::pair<std::string, std::string>> &jobs,
                       const std::string &skimConfig,
                       std::size_t &nextJobIndex,
                       std::size_t &written,
                       std::size_t batchLimit)
{
    while (nextJobIndex < jobs.size())
    {
        const auto &job = jobs[nextJobIndex];
        const auto &inputFile = job.first;
        const auto &outFile = job.second;
        ofs << inputFile << " " << outFile << " " << skimConfig << "\n";
        ++nextJobIndex;
        ++written;
        if (batchLimit > 0 && written >= batchLimit)
            return true;
    }
    return true;
}
} // namespace

int main(int argc, char *argv[])
{
    std::string skimConfig;
    std::string jobOutBase;
    std::size_t batchLimit = 1500;

    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        if (arg == "--skimConfig" && i + 1 < argc)
        {
            skimConfig = argv[++i];
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
            rdfWS_utility::messageWARN("skimSamples_slurm",
                                       "Unknown or incomplete argument: " + arg);
            return 1;
        }
    }

    if (skimConfig.empty() || jobOutBase.empty())
    {
        rdfWS_utility::messageERROR(
            "skimSamples_slurm",
            "Input template wrong. Example: \"skimSamples_slurm --skimConfig <skim json> --jobOut <job base> [--batchLimit N]\"");
        return 1;
    }

    auto jsonData = rdfWS_utility::readJson("skimSamples_slurm", skimConfig);
    rdfWS_utility::JsonObject configFile(jsonData, "Skim JO config");

    if (configFile.contains("jobType"))
    {
        std::string jobType = configFile.at("jobType").get<std::string>();
        if (jobType != "skim")
        {
            rdfWS_utility::messageERROR("skimSamples_slurm", "The jobType of your config is not skimming!");
            return 1;
        }
    }

    int isData = configFile.at("isData").get<int>();
    std::string era = configFile.at("era").get<std::string>();
    std::string outDir = configFile.at("outDir").get<std::string>();

    int maxFiles = 0;
    if (configFile.contains("maxFiles"))
    {
        maxFiles = configFile.at("maxFiles").get<int>();
        if (maxFiles < 0)
            maxFiles = 0;
    }

    std::vector<std::string> channels = configFile.at("datasets").get<std::vector<std::string>>();
    std::string sampleConfigPath = configFile.at("sampleConfig").get<std::string>();
    if (sampleConfigPath.empty())
        sampleConfigPath = "json/samples/" + era + ".json";

    SampleControl samples(sampleConfigPath);

    const std::filesystem::path outPath(jobOutBase);
    if (!outPath.parent_path().empty() && !std::filesystem::exists(outPath.parent_path()))
    {
        std::filesystem::create_directories(outPath.parent_path());
    }

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
            rdfWS_utility::messageERROR("skimSamples_slurm", "Failed to open job file: " + jobOut);
            return false;
        }
        return true;
    };
    if (!openBatch())
        return 1;

    for (const auto &channel : channels)
    {
        auto filePaths = samples.getFiles(channel);
        if (filePaths.empty())
            continue;

        std::vector<std::pair<std::string, std::string>> jobs;
        jobs.reserve(filePaths.size());
        for (const auto &filePath : filePaths)
        {
            std::string outputPath = buildOutputPath(outDir, isData, era, channel, filePath);
            if (std::filesystem::exists(outputPath))
            {
                TFile fcheck(outputPath.c_str(), "READ");
                if (fcheck.IsZombie())
                {
                    rdfWS_utility::messageWARN("skimSamples_slurm", "Output is zombie, will reskim: " + outputPath);
                    fcheck.Close();
                }
                else if (!fcheck.Get("genWeightSum"))
                {
                    rdfWS_utility::messageWARN("skimSamples_slurm", "Output missing genWeightSum, will reskim: " + outputPath);
                    fcheck.Close();
                }
                else
                {
                    rdfWS_utility::messageINFO("skimSamples_slurm", outputPath + " already exists and is valid, skip");
                    continue;
                }
            }
            jobs.emplace_back(filePath, outputPath);
        }

        if (maxFiles > 0 && jobs.size() > static_cast<std::size_t>(maxFiles))
        {
            jobs.erase(jobs.begin() + maxFiles, jobs.end());
        }

        std::size_t nextJobIndex = 0;
        while (nextJobIndex < jobs.size())
        {
            if (!appendJobs(ofs, jobs, skimConfig, nextJobIndex, writtenInBatch, batchLimit))
                return 1;
            if (nextJobIndex >= jobs.size() || batchLimit == 0 || writtenInBatch < batchLimit)
                break;
            ++batchIndex;
            writtenInBatch = 0;
            if (!openBatch())
                return 1;
        }
    }

    return 0;
}
