#include "SampleControl.h"
#include "Utility.h"

#include "external/json.hpp"

#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <fstream>

SampleControl::SampleControl(const std::string &jsonFileName)
{
    // reading data json file path
    // e.g. samples/Run2023D.json
    auto jsonData = rdfWS_utility::readJson("SampleControl", jsonFileName);

    // "dir" and "file" entry are necessary
    if (!jsonData.contains("dir") || !jsonData.contains("file"))
    {
        rdfWS_utility::messageERROR("SampleControl", "JSON file " + jsonFileName + " lack required keys: 'dir' or 'file'");
    }

    // extract dataset infor
    std::map<std::string, std::string> dirMap;
    std::map<std::string, std::vector<std::string>> samplesMap;
    try
    {
        dirMap = jsonData["dir"].get<std::map<std::string, std::string>>();
        samplesMap = jsonData["file"].get<std::map<std::string, std::vector<std::string>>>();
    }
    catch (const nlohmann::json::type_error &ex)
    {
        rdfWS_utility::messageERROR("SampleControl", "JSON file " + jsonFileName + " structure error: " + std::string(ex.what()));
    }

    // check dataset and record
    for (const auto &[key, dir] : dirMap)
    {
        if (!std::filesystem::exists(dir))
        {
            rdfWS_utility::messageWARN("SampleControl", "Directory " + dir + " does not exist.");
            continue;
        }
        if (samplesMap.find(key) == samplesMap.end())
        {
            rdfWS_utility::messageWARN("SampleControl", "Missing directory " + key + " for file block.");
            continue;
        }
        this->_types.push_back(key);
        this->_filePaths.emplace(key, std::vector<std::string>{});
        this->_missingFiles.emplace(key, std::vector<std::string>{});
        // record valid files
        const auto &baseDir = dirMap[key];
        for (const auto &file : samplesMap[key])
        {
            std::string fullPath = baseDir + "/" + file;
            if (!std::filesystem::exists(fullPath))
            {
                this->_missingFiles[key].push_back(fullPath);
                rdfWS_utility::messageWARN("SampleControl", "File " + file + " does not exist in directory " + dir + ".");
            }
            else
            {
                this->_filePaths[key].push_back(fullPath);
            }
        }
        this->_isValid[key] = _missingFiles[key].empty() ? 1 : 0;
    }
}

std::vector<std::string> SampleControl::getFiles(std::string type)
{
    if (_isValid[type] == 0)
    {
        rdfWS_utility::messageWARN("SampleControl", "Files missing: ");
        auto missingFiles = this->_missingFiles[type];
        for (auto file : missingFiles)
        {
            rdfWS_utility::messageWARN("SampleControl", "\t"+file);
        }
    }

    if (this->_filePaths.find(type) == this->_filePaths.end())
    {
        rdfWS_utility::messageWARN("SampleControl", "Requested type '" + type + "' does not exist in loaded samples.");
    }

    return this->_filePaths[type];
}

std::vector<std::string> SampleControl::getMergedFiles(const std::vector<std::string> &types)
{
    std::vector<std::string> mergedFiles;

    for (const auto &type : types)
    {
        if (this->_filePaths.find(type) == this->_filePaths.end())
        {
            rdfWS_utility::messageWARN("SampleControl", "Warning: Type " + type + " not found. Will skip");
            continue;
        }

        // Append files for the current type
        mergedFiles.insert(mergedFiles.end(), this->_filePaths[type].begin(), this->_filePaths[type].end());
    }

    if (mergedFiles.empty())
    {
        rdfWS_utility::messageWARN("SampleControl", "No valid files found for the provided types.");
    }

    return mergedFiles;
}
