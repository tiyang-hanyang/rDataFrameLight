#ifndef SAMPLE_CONTROL_H
#define SAMPLE_CONTROL_H

#include <string>
#include <vector>
#include <map>

class SampleControl
{
    private:
        std::vector<std::string> _types;
        std::map<std::string, std::vector<std::string>> _filePaths;
        std::map<std::string, int> _isValid; // 0: files incomplete, 1: good
        std::map<std::string, std::vector<std::string>> _missingFiles;
    public:
        explicit SampleControl(const std::string& jsonFileName);
        std::vector<std::string> getFiles(std::string type);
        std::vector<std::string> getMergedFiles(const std::vector<std::string>& types);
};


#endif 