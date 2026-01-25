#ifndef RDATAFRAME_LIGHT_UTILITY_H
#define RDATAFRAME_LIGHT_UTILITY_H

#include "external/json.hpp"

#include <string>

/// @brief A shared namespace, providing common utilities needed by various control modules
namespace rdfWS_utility
{
    nlohmann::json readJson(const std::string &prog, const std::string &configPath);

    void messageINFO(const std::string &prog, const std::string &content);
    void messageWARN(const std::string &prog, const std::string &content);
    [[noreturn]] void messageERROR(const std::string &prog, const std::string &content);

    void creatingFolder(const std::string &prog, const std::string &folderDir);

    template <typename T>
    T parseJson(nlohmann::json json, std::string jsonName, std::string term, std::string progName)
    {
        try
        {
            T readObject = json[term];
            return readObject;
        }
        catch (const nlohmann::json::exception &e)
        {
            messageERROR(progName, "Reading term " + term + " from the json " + jsonName + " failed: " + e.what());
        }
    };

    class JsonObject
    {
    private:
        nlohmann::json _json;
        std::string _name;

    public:
        JsonObject(nlohmann::json jsonConfig, const std::string &jsonName);

        bool contains(const std::string& key);
        JsonObject at(const std::string &item);

        template <typename T>
        T get() const
        {
            try
            {
                return this->_json.get<T>();
            }
            catch (const nlohmann::json::exception &e)
            {
                messageERROR(this->_name, "JSON " + this->_name + " error: " + e.what());
            }
        }

        template <typename T>
        operator T()
        {
            try
            {
                return this->_json.get<T>();
            }
            catch (const nlohmann::json::exception &e)
            {
                messageERROR(this->_name, "JSON " + this->_name + " error: " + e.what());
            }
        }
    };

}

#endif