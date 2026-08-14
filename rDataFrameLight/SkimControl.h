#ifndef SKIM_CONTROL_H
#define SKIM_CONTROL_H

#include "CutControl.h"
#include "SampleControl.h"

#include "external/json.hpp"

#include <vector>
#include <map>
#include <string>
#include <set>
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

    // enable to store in separated or merging samples
    int _ifMerging;

    // the very top outDir folder, if not exist, create it
    std::string _outDir;

    // all the channels processed in this run
    std::vector<std::string> _channels;
    // std::map<std::string, int> _isData;
    int _isData;
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

    // in case distinguish preliminary selection and further cut
    int _isPreliminary;
    // to deal with the stack up samples.
    int _isStack;

    // limit number of files per channel (0 means no limit)
    int _maxFilesPerChannel;

    // Optional syst-aware loose skim filter, intended to keep events that can
    // pass the 4-jet selection under the maximum JME up variation.
    int _useSystAwareFourJet = 0;
    std::string _systAwareSelectionType = "fourJetUpEnvelope";
    int _systAwareNJet = 4;
    float _systAwareJetPtThreshold = 30.0;
    float _systAwareBTagThreshold = 0.1272;
    std::string _systAwareBTagBranch = "Jet_btagUParTAK4B";
    std::string _systAwareGoodJetExpr = "";
    int _systAwareRequireJVMEnvelope = 1;
    int _skipNominalFourJetSteps = 1;
    std::map<std::string, std::vector<std::string>> _skimShiftSyst;
    std::set<std::string> _skimSystSkipSteps;

public:
    SkimControl() = default;
    explicit SkimControl(nlohmann::json configFile);
    explicit SkimControl(const std::string& configPath);
    ~SkimControl() = default;

    void readConfig(nlohmann::json configFile);
    void readConfig(const std::string& configPath);

    // turn off & turn on channels to run
    void turnOn(const std::string& channels);
    void turnOff(const std::string& channels);

    // operations inside the run
    double _getTotalGenWeight(std::vector<std::string> fileLists);
    ROOT::RDF::RNode _preliminaryDeco(ROOT::RDF::RNode rndDS, const std::string& channel, double totalGenWeight);
    ROOT::RDF::RNode _applySystAwareFourJet(ROOT::RDF::RNode rndDS);
    std::vector<std::string> _getBranchArray(ROOT::RDF::RNode rndDS, int isPreliminary);

    // better to split run one by one file, to avoid failure in the middle
    void run();
};

#endif
