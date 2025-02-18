#include "Utility.h"

#include <iostream>
#include <filesystem>
#include <fstream>
#include <stdexcept>

/// @brief Function to safely read the json for different control modules
/// @param prog The name of the module, enabling easy debugging.
/// @param jsonPath The path of the target json
/// @return The json format data structure, details checking and extraction is left to specific modules.
nlohmann::json rdfWS_utility::readJson(const std::string &prog, const std::string &configPath)
{
    // if not absolute path, need to make sure json is in reference of the source/ folder
    std::string properPath = configPath;
    if (configPath[0] != '/')
    {
        std::filesystem::path sourceFile = __FILE__;
        properPath = sourceFile.parent_path().parent_path().string() + "/" + configPath;
    }

    messageINFO(prog, "Loading json from: " + properPath);
    try
    {
        // check existence of the json file
        if (!std::filesystem::exists(properPath))
        {
            messageERROR(prog, "JSON file " + properPath + " not exists");
        }
    }
    catch (const std::filesystem::filesystem_error &ex)
    {
        // handle possibly catch exception due to authorization issue or disk problems
        messageERROR(prog, "Filesystem error: " + std::string(ex.what()));
    }

    // load json
    nlohmann::json jsonData;
    try
    {
        std::ifstream inputFile(properPath);
        inputFile >> jsonData;
    }
    catch (const nlohmann::json::parse_error &ex)
    {
        messageERROR(prog, "JSON file " + properPath + " parse error: " + std::string(ex.what()));
    }

    return jsonData;
}

/// @brief INFO level message.
/// @param prog The name of the module, enabling easy debugging.
/// @param content
void rdfWS_utility::messageINFO(const std::string &prog, const std::string &content)
{
    std::cout << "[" << prog << "] INFO: " << content << std::endl;
}

/// @brief WARNING level message.
/// @param prog The name of the module, enabling easy debugging.
/// @param content
void rdfWS_utility::messageWARN(const std::string &prog, const std::string &content)
{
    std::cerr << "[" << prog << "] WARNING: " << content << std::endl;
}

/// @brief ERROR level message, will throw.
/// @param prog The name of the module, enabling easy debugging.
/// @param content
[[noreturn]] void rdfWS_utility::messageERROR(const std::string &prog, const std::string &content)
{
    throw std::runtime_error("[" + prog + "] ERROR: " + content);
}

/// @brief Function to safely create a folder.
/// @param prog The name of the module, enabling easy debugging.
/// @param folderDir The directory of the folder.
void rdfWS_utility::creatingFolder(const std::string &prog, const std::string &folderDir)
{
    // mkdir for output
    if (!std::filesystem::exists(folderDir))
    {
        if (std::filesystem::create_directories(folderDir))
        {
            messageINFO(prog, "Dir " + folderDir + " created.");
        }
        else
        {
            messageERROR(prog, "Dir " + folderDir + " created fail! Please check your authorization or disk quota!");
        }
    }
    else
    {
        messageINFO(prog, "Dir " + folderDir + " already exists.");
    }
}