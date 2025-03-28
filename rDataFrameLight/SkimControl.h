#ifndef SKIM_CONTROL_H
#define SKIM_CONTROL_H

#include "CutControl.h"
#include "SampleControl.h"

#include "external/json.hpp"

#include <vector>
#include <map>
#include <string>
#include <functional>
#include <optional>

#include "ROOT/RDataFrame.hxx"

// it is hard to have a general reweighting tool
// thus I put the reweighting just inside this SkimControl
class SkimControl
{
private:
    // the name of this skimming job
    std::string _skimName;
    // year and era to determine the golden json and MC campaign 
    std::string _run;
    std::string _year;
    std::string _era;

    // the very top outDir folder, if not exist, create it
    std::string _outDir;

    // all the channels processed in this run
    std::vector<std::string> _channels;
    std::map<std::string, int> _isData;
    // for mc, the XS value of this samples need to be computed
    std::string _mcWeight;
    std::map<std::string, float> _XSvalues;
    std::optional<SampleControl> _samples;
    // allow allow manually turn off
    std::map<std::string, int> _isOn;

    // branch to keep
    std::vector<std::string> _branchList;

    // Encode the lambda for applying the golden json for reusage
    std::optional<std::function<ROOT::RDF::RNode(ROOT::RDF::RNode)>> _goldenJsonLambda;

    // The main filter of the skim
    CutControl _skimCut;

    // member function to apply the golden json
    void _createGoldenJsonFunc();
    // ROOT::RDF::RNode applyGoldenJson(ROOT::RDF::RNode origData);

    // for gracefully exit
    static void signalHandler(int signum);
    static SkimControl* instance;
    std::atomic<bool> stop_requested{false};

public:
    SkimControl() = default;
    explicit SkimControl(nlohmann::json configFile);
    explicit SkimControl(const std::string& configPath);
    ~SkimControl() = default;
    // no need to copy skim jobs I suppose

    void readConfig(nlohmann::json configFile);
    void readConfig(const std::string& configPath);

    // turn off & turn on channels to run
    void turnOn(const std::string& channels);
    void turnOff(const std::string& channels);

    // better to split run one by one file, to avoid failure in the middle
    void run();
};

#endif