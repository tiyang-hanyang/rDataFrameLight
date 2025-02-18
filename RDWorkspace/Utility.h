#ifndef RDATAFRAME_LIGHT_UTILITY_H
#define RDATAFRAME_LIGHT_UTILITY_H

#include "external/json.hpp"

#include <string>

/// @brief A shared namespace, providing common utilities needed by various control modules
namespace rdfWS_utility
{
    nlohmann::json readJson(const std::string& prog, const std::string& configPath);

    void messageINFO(const std::string& prog, const std::string& content);
    void messageWARN(const std::string& prog, const std::string& content);
    [[noreturn]] void messageERROR(const std::string& prog, const std::string& content);

    void creatingFolder(const std::string& prog, const std::string& folderDir);
}

#endif