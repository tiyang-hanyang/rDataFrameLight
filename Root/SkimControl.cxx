#include "SkimControl.h"
#include "Utility.h"

#include <filesystem>

namespace
{
std::map<std::string, std::vector<std::string>> loadSystAliases()
{
    rdfWS_utility::JsonObject aliasConfig(
        rdfWS_utility::readJson("SkimControl", "json/general_config/syst_NP.json"),
        "Systematic Alias Config"
    );
    return aliasConfig.get<std::map<std::string, std::vector<std::string>>>();
}

std::map<std::string, std::vector<std::string>> expandShiftSystAliases(
    const nlohmann::json &rawConfig)
{
    if (!rawConfig.is_object())
        rdfWS_utility::messageERROR("SkimControl", "shiftSyst must be a map.");

    const auto aliases = loadSystAliases();
    std::map<std::string, std::vector<std::string>> expanded;
    for (auto it = rawConfig.begin(); it != rawConfig.end(); ++it)
    {
        if (!it.value().is_array())
            rdfWS_utility::messageERROR("SkimControl", "shiftSyst entry " + it.key() + " must be a list of target branches.");
        const auto targets = it.value().get<std::vector<std::string>>();
        auto aliasIter = aliases.find(it.key());
        if (aliasIter != aliases.end())
        {
            for (const auto &name : aliasIter->second)
                expanded.emplace(name, targets);
        }
        else
        {
            expanded.emplace(it.key(), targets);
        }
    }
    return expanded;
}

}
#include <fstream>
#include <stdexcept>
#include <memory>
#include <regex>
#include <csignal>
#include <unistd.h>
#include <atomic>
#include <unordered_map>
#include <algorithm>
#include <cctype>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TInterpreter.h"
#include "ROOT/RSnapshotOptions.hxx"
#include "ROOT/RDFHelpers.hxx"

namespace
{
    void DeclareSkimSystHelpersToCling()
    {
        static bool done = false;
        if (done)
            return;
        gInterpreter->Declare(R"(
            #include "ROOT/RVec.hxx"
            #include <algorithm>
            #include <cmath>
            #include <initializer_list>
            #include <stdexcept>

            ROOT::VecOps::RVec<float> skimQuadratureSum(std::initializer_list<ROOT::VecOps::RVec<float>> sources)
            {
                ROOT::VecOps::RVec<float> total;
                bool initialized = false;
                for (const auto& source : sources)
                {
                    if (!initialized)
                    {
                        total = ROOT::VecOps::RVec<float>(source.size(), 0.0f);
                        initialized = true;
                    }
                    if (source.size() != total.size())
                    {
                        throw std::runtime_error("skimQuadratureSum got RVec inputs with different sizes.");
                    }
                    for (size_t i = 0; i < source.size(); ++i)
                    {
                        total[i] += source[i] * source[i];
                    }
                }
                for (size_t i = 0; i < total.size(); ++i)
                {
                    total[i] = std::sqrt(total[i]);
                }
                return total;
            }

            ROOT::VecOps::RVec<float> skimJMEUpEnvelopePt(
                const ROOT::VecOps::RVec<float>& Jet_pt_JEC,
                const ROOT::VecOps::RVec<float>& CMS_scale_j_Total,
                const ROOT::VecOps::RVec<float>& Jet_pt_JEScorr,
                const ROOT::VecOps::RVec<float>& Jet_JER_corr_up,
                const ROOT::VecOps::RVec<float>& Jet_JER_corr_down)
            {
                int nJet = Jet_pt_JEC.size();
                ROOT::VecOps::RVec<float> shifted(nJet);
                for (int i = 0; i < nJet; ++i)
                {
                    const float jesAbs = Jet_pt_JEC[i] * CMS_scale_j_Total[i];
                    const float jerUpPt = Jet_pt_JEScorr[i] * Jet_JER_corr_up[i];
                    const float jerDownPt = Jet_pt_JEScorr[i] * Jet_JER_corr_down[i];
                    const float jerMaxPt = std::max(jerUpPt, jerDownPt);
                    const float jerAbs = std::max(0.0f, jerMaxPt - Jet_pt_JEC[i]);
                    shifted[i] = Jet_pt_JEC[i] + std::sqrt(jesAbs * jesAbs + jerAbs * jerAbs);
                }
                return shifted;
            }

            ROOT::VecOps::RVec<float> skimMinDistanceFromMuon(
                const ROOT::VecOps::RVec<int>& isGoodJet,
                const ROOT::VecOps::RVec<float>& eta,
                const ROOT::VecOps::RVec<float>& phi,
                const int& leadingMuonIdx,
                const int& subleadingMuonIdx,
                const ROOT::VecOps::RVec<float>& Muon_eta,
                const ROOT::VecOps::RVec<float>& Muon_phi)
            {
                auto jetSize = isGoodJet.size();
                ROOT::VecOps::RVec<float> minDR(jetSize);
                auto leadingMuonEta = Muon_eta[leadingMuonIdx];
                auto leadingMuonPhi = Muon_phi[leadingMuonIdx];
                auto subleadingMuonEta = Muon_eta[subleadingMuonIdx];
                auto subleadingMuonPhi = Muon_phi[subleadingMuonIdx];
                for (auto i=0; i < jetSize; i++)
                {
                    if (!isGoodJet[i])
                    {
                        minDR[i]=0.0;
                        continue;
                    }
                    float dRleading = ROOT::VecOps::DeltaR(leadingMuonEta, eta[i], leadingMuonPhi, phi[i]);
                    float dRsubleading = ROOT::VecOps::DeltaR(subleadingMuonEta, eta[i], subleadingMuonPhi, phi[i]);
                    minDR[i] = std::min(dRleading, dRsubleading);
                }
                return minDR;
            }
        )");
        done = true;
    }
}

////////////////////////////////////////////////// Setup configs

SkimControl::SkimControl(nlohmann::json configFile)
{
    this->readConfig(configFile);
}

SkimControl::SkimControl(const std::string &configPath)
{
    this->readConfig(configPath);
}

void SkimControl::readConfig(nlohmann::json origConfigFile)
{
    rdfWS_utility::JsonObject configFile(origConfigFile, "Skim JO config");

    // job type check, better to have to avoid confusion
    if (configFile.contains("jobType"))
    {
        std::string jobType = configFile.at("jobType").get<std::string>();
        if (jobType != "skim") 
        {
            rdfWS_utility::messageERROR("SkimControl.cxx", "The jobType of your config is not skimming! Please check again. Running ceases...");
            exit(1);
        }
    }

    this->_ifMerging = configFile.at("merge").get<int>();
    this->_year = configFile.at("year").get<std::string>();
    this->_era = configFile.at("era").get<std::string>();
    this->_skimName = this->_era + "_skimmmed";
    if (configFile.contains("name"))
    {
       this->_skimName = configFile.at("name").get<std::string>();
    }

    // by default do preliminary selection
    // mc sum genWeight from Events tree
    // data need test the golden json
    int preliminary(1);
    if (configFile.contains("preliminary"))
    {
        if (configFile.at("preliminary").get<int>() == 0) preliminary=0;
    }
    this->_isPreliminary = preliminary;

    // set maximum files for skimming
    int maxFiles(0);
    if (configFile.contains("maxFiles"))
    {
        maxFiles = configFile.at("maxFiles").get<int>();
        if (maxFiles < 0) maxFiles = 0;
    }
    this->_maxFilesPerChannel = maxFiles;

    if (origConfigFile.contains("skimSystSelection"))
    {
        auto systSelection = origConfigFile.at("skimSystSelection");
        std::string selectionType = systSelection.value("type", "");
        if (selectionType == "fourJetUpEnvelope" || selectionType == "bJetUpEnvelope")
        {
            this->_useSystAwareFourJet = 1;
            this->_systAwareSelectionType = selectionType;
            this->_systAwareNJet = systSelection.value("nJet", selectionType == "bJetUpEnvelope" ? 3 : 4);
            this->_systAwareJetPtThreshold = systSelection.value("jetPtThreshold", 30.0);
            this->_systAwareBTagThreshold = systSelection.value("btagThreshold", 0.1272);
            this->_systAwareBTagBranch = systSelection.value("btagBranch", "Jet_btagUParTAK4B");
            this->_systAwareGoodJetExpr = systSelection.value("goodJetExpr", "");
            this->_systAwareRequireJVMEnvelope = systSelection.value("requireJVMEnvelope", selectionType == "fourJetUpEnvelope" ? 1 : 0);
            this->_skipNominalFourJetSteps = systSelection.value("skipNominalFourJet", 1);
            if (systSelection.contains("skipSteps"))
            {
                auto skipSteps = systSelection.at("skipSteps").get<std::vector<std::string>>();
                this->_skimSystSkipSteps.insert(skipSteps.begin(), skipSteps.end());
            }
            else
            {
                this->_skimSystSkipSteps.insert("pass JVM necessary");
                if (selectionType == "bJetUpEnvelope")
                {
                    this->_skimSystSkipSteps.insert("BJetCond");
                    this->_skimSystSkipSteps.insert("nBJet");
                    this->_skimSystSkipSteps.insert("at least three bjet");
                }
                else
                {
                    this->_skimSystSkipSteps.insert("nGoodJet");
                    this->_skimSystSkipSteps.insert("at least four good jets");
                }
            }
            if (this->_systAwareRequireJVMEnvelope)
                this->_skimSystSkipSteps.insert("pass JVM necessary");

            if (systSelection.contains("shiftSyst"))
            {
                this->_skimShiftSyst = expandShiftSystAliases(systSelection.at("shiftSyst"));
            }
            else if (origConfigFile.contains("shiftSyst"))
            {
                this->_skimShiftSyst = expandShiftSystAliases(origConfigFile.at("shiftSyst"));
            }
            else
            {
                rdfWS_utility::messageERROR("SkimControl", "skimSystSelection " + selectionType + " requires shiftSyst either inside skimSystSelection or at top level.");
            }
        }
        else if (selectionType != "")
        {
            rdfWS_utility::messageERROR("SkimControl", "Unknown skimSystSelection type: " + selectionType);
        }
    }

    this->_outDir = configFile.at("outDir").get<std::string>();
    rdfWS_utility::creatingFolder("SkimControl", this->_outDir);

    // channel info
    this->_channels = configFile.at("datasets").get<std::vector<std::string>>();
    // by default, turn on all the channel skim
    for (const auto &channel : this->_channels)
    {
        this->_isOn.emplace(channel, 1);
    }
    // this->_isData = configFile.at("isData").get<std::map<std::string, int>>();
    this->_isData = configFile.at("isData").get<int>();

    // files
    std::string sampleConfigPath = configFile.at("sampleConfig");
    if (sampleConfigPath == "")
        this->_samples.emplace("json/samples/" + this->_era + ".json");
    else
        this->_samples.emplace(sampleConfigPath);

    // prepare all the cuts
    std::vector<std::string> cutConfigList = configFile.at("cutConfig");
    if (cutConfigList.size() > 0)
    {
        this->_skimCut = CutControl(cutConfigList[0]);
        for (int i = 1; i < cutConfigList.size(); i++)
        {
            this->_skimCut = this->_skimCut + CutControl(cutConfigList[i]);
        }
    }

    // branch list
    std::string branchPath = configFile.at("branchConfig");
    rdfWS_utility::JsonObject branchJson(rdfWS_utility::readJson("SkimControl", branchPath), "Branch Config");
    this->_branchList = branchJson.get<std::vector<std::string>>();
}

void SkimControl::readConfig(const std::string &configPath)
{
    auto jsonData = rdfWS_utility::readJson("SkimControl", configPath);
    this->readConfig(jsonData);
}

////////////////////////////////////////////////// Setup Data GoldenJson Lambda
void SkimControl::_createGoldenJsonFunc()
{
    rdfWS_utility::messageINFO("SkimControl", "Creating the golden json lambda function.");
    // check the golden json re-assignment
    if (this->_goldenJsonLambda)
    {
        rdfWS_utility::messageWARN("SkimControl", "The _goldenJsonLambda is already initialized. It is re-initialized again. Please make sure this is really what desired.");
    }

    // parsing the goldenJson, cast from string to int
    std::string goldenJsonPath = "json/goldenJson/" + this->_year + "/goldenJson.json";
    // rdfWS_utility::messageINFO("SkimControl", "Loading golden json from: " + goldenJsonPath);
    auto goldenJson = rdfWS_utility::readJson("SkimControl", goldenJsonPath);
    std::map<std::string, std::vector<std::pair<unsigned int, unsigned int>>> goldenJsonRaw = goldenJson;
    std::unordered_map<unsigned int, std::vector<std::pair<unsigned int,unsigned int>>> goldenJsonList;
    goldenJsonList.reserve(goldenJsonRaw.size());

    auto is_digits = [](const std::string& s) {
        return !s.empty() && std::all_of(s.begin(), s.end(), [](unsigned char c) { return std::isdigit(c) != 0; });
    };

    for (const auto& [k,v] : goldenJsonRaw) 
    {
        if (!is_digits(k)) {
            rdfWS_utility::messageWARN("SkimControl", "Invalid run key: " + k + ", skip.");
            continue;
        }
        unsigned int run = std::stoul(k);
        goldenJsonList.emplace(run, v);
    }

    // golden json function definition
    this->_goldenJsonLambda =
        [goldenJsonList](ROOT::RDF::RNode origData)
    { return origData.Filter(
        [goldenJsonList](unsigned int run, unsigned int luminosityBlock)
        {
            if (goldenJsonList.find(run) == goldenJsonList.end())
                return false;
            auto goldenLumiBlk = goldenJsonList.at(run);
            int isGood(0);
            for (const auto &blkPair : goldenLumiBlk)
            {
                if (luminosityBlock <= blkPair.second && luminosityBlock >= blkPair.first)
                {
                    isGood = 1;
                    break;
                }
            }
            if (!isGood)
                return false;
            return true;
        },
        {"run", "luminosityBlock"}); };
}

void SkimControl::turnOn(const std::string &channels)
{
    bool found = false;

    std::regex pattern(channels);
    try
    {
        std::regex pattern(channels);

        for (auto &[key, value] : this->_isOn)
        {
            if (std::regex_match(key, pattern))
            {
                value = 1;
                found = true;
                rdfWS_utility::messageINFO("SkimControl", "Turned on channel: " + key);
            }
        }
    }
    catch (const std::regex_error &)
    {
        rdfWS_utility::messageWARN("SkimControl", "Invalid channel pattern: " + channels);
        return;
    }

    if (!found)
    {
        rdfWS_utility::messageWARN("SkimControl", "No channels matched: " + channels);
    }
}

void SkimControl::turnOff(const std::string &channels)
{
    bool found = false;

    std::regex pattern(channels);
    try
    {
        std::regex pattern(channels);

        for (auto &[key, value] : this->_isOn)
        {
            if (std::regex_match(key, pattern))
            {
                value = 0;
                found = true;
                rdfWS_utility::messageINFO("SkimControl", "Turned off channel: " + key);
            }
        }
    }
    catch (const std::regex_error &)
    {
        rdfWS_utility::messageWARN("SkimControl", "Invalid channel pattern: " + channels);
        return;
    }

    if (!found)
    {
        rdfWS_utility::messageWARN("SkimControl", "No channels matched: " + channels);
    }
}

////////////////////////////////////////////////// Enable exiting gracefully

SkimControl *SkimControl::instance = nullptr;

void SkimControl::signalHandler(int signum)
{
    if (instance)
    {
        instance->stop_requested = true;
    }
}

////////////////////////////////////////////////// Skimming methods
// function for getting the total weight if exist
double SkimControl::_getTotalGenWeight(std::vector<std::string> fileLists)
{
    double totalGenWeight(0.0);
    for (const auto& fPath: fileLists)
    {
        std::unique_ptr<TFile> fTemp{ TFile::Open(fPath.c_str(), "READ") };
        if (!fTemp || fTemp->IsZombie()) 
        {
            rdfWS_utility::messageWARN("SkimControl", fPath+" not exist!");
            continue;
        }

        auto sumWeightHist = dynamic_cast<TH1D*>(fTemp->Get("genWeightSum"));
        if (!sumWeightHist)
        {
            rdfWS_utility::messageWARN("SkimControl", fPath+" does not have genWeightSum! Try manually add.");
            if (!fTemp->GetListOfKeys()->FindObject("Events")) continue;
            auto tTemp = dynamic_cast<TTree*>(fTemp->Get("Events"));
            if (!tTemp) 
            {
                rdfWS_utility::messageWARN("SkimControl", fPath+" Events is not a TTree, skip...");
                continue;
            }

            ROOT::RDataFrame dfTemp(*tTemp);
            try {
                totalGenWeight += dfTemp.Sum("genWeight").GetValue();
            } catch (const std::exception& e) {
                rdfWS_utility::messageWARN("SkimControl", fPath+" cannot add genWeight, skip...");
                continue;
            }
        }
        else
            totalGenWeight += sumWeightHist->GetBinContent(1);
    }
    rdfWS_utility::messageINFO("SkimControl", "Total genWeight: " + std::to_string(totalGenWeight) );
    return totalGenWeight;
}

ROOT::RDF::RNode SkimControl::_preliminaryDeco(ROOT::RDF::RNode rndDS, const std::string& channel, double totalGenWeight) 
{
    // apply golden json
    if (!(this->_goldenJsonLambda))
    {
        this->_createGoldenJsonFunc();
    }
    rndDS = this->_goldenJsonLambda.value()(rndDS);

    return rndDS;
}

ROOT::RDF::RNode SkimControl::_applySystAwareFourJet(ROOT::RDF::RNode rndDS)
{
    if (!this->_useSystAwareFourJet)
        return rndDS;

    DeclareSkimSystHelpersToCling();

    if (this->_systAwareRequireJVMEnvelope)
    {
        std::string jvmEnvelopeExpr = "(JVMweight > 0.0";
        auto availableColumns = rndDS.GetColumnNames();
        for (const auto &column : availableColumns)
        {
            if (column.rfind("JVMweight_", 0) == 0)
                jvmEnvelopeExpr += " || " + column + " > 0.0";
        }
        jvmEnvelopeExpr += ")";
        rndDS = rndDS.Define("JVMweight_skimEnvelope", jvmEnvelopeExpr);
    }

    std::vector<std::string> jetPtSources;
    for (const auto &[systName, targets] : this->_skimShiftSyst)
    {
        if (std::find(targets.begin(), targets.end(), "Jet_pt_JEC") != targets.end() && rndDS.HasColumn(systName))
            jetPtSources.push_back(systName);
    }
    if (jetPtSources.empty())
    {
        rdfWS_utility::messageERROR("SkimControl", "syst-aware skim needs at least one shiftSyst source targeting Jet_pt_JEC.");
    }

    std::string totalExpr = "skimQuadratureSum({";
    for (size_t i = 0; i < jetPtSources.size(); ++i)
    {
        if (i > 0)
            totalExpr += ", ";
        totalExpr += jetPtSources[i];
    }
    totalExpr += "})";

    rndDS = rndDS.Define("CMS_scale_j_skimTotalUpEnvelope", totalExpr);
    if (rndDS.HasColumn("Jet_pt_JEScorr") && rndDS.HasColumn("Jet_JER_corr_up") && rndDS.HasColumn("Jet_JER_corr_down"))
    {
        rndDS = rndDS.Define("Jet_pt_JEC_skimJMEUpEnvelope", "skimJMEUpEnvelopePt(Jet_pt_JEC, CMS_scale_j_skimTotalUpEnvelope, Jet_pt_JEScorr, Jet_JER_corr_up, Jet_JER_corr_down)");
    }
    else
    {
        rdfWS_utility::messageWARN("SkimControl", "JER columns not available; syst-aware fourJet skim uses JES-only up envelope.");
        rndDS = rndDS.Define("Jet_pt_JEC_skimJMEUpEnvelope", "Jet_pt_JEC * (1.0f + CMS_scale_j_skimTotalUpEnvelope)");
    }

    if (this->_skipNominalFourJetSteps)
        return rndDS;

    if (this->_systAwareGoodJetExpr.empty())
    {
        if (rndDS.HasColumn("Jet_tightNoPt") && rndDS.HasColumn("Jet_drFromMuon"))
        {
            rndDS = rndDS.Define(
                "GoodJetCond_skimJMEUpEnvelope",
                "(Jet_pt_JEC_skimJMEUpEnvelope > " + std::to_string(this->_systAwareJetPtThreshold) + ") && Jet_tightNoPt && (Jet_drFromMuon > 0.4)"
            );
        }
        else
        {
            rndDS = rndDS.Define(
                "Jet_mediumPtTight_skimJMEUpEnvelope",
                "(Jet_pt_JEC_skimJMEUpEnvelope > " + std::to_string(this->_systAwareJetPtThreshold) + ") && (abs(Jet_eta)<2.5) && (Jet_rawFactor<0.9) && JetIdTight"
            );
            rndDS = rndDS.Define(
                "Jet_drFromMuon_skimJMEUpEnvelope",
                "skimMinDistanceFromMuon(Jet_mediumPtTight_skimJMEUpEnvelope, Jet_eta, Jet_phi, leadingMuonIdx, subleadingMuonIdx, Muon_eta, Muon_phi)"
            );
            rndDS = rndDS.Define("GoodJetCond_skimJMEUpEnvelope", "(Jet_drFromMuon_skimJMEUpEnvelope>0.4)");
        }
    }
    else
    {
        rndDS = rndDS.Define("GoodJetCond_skimJMEUpEnvelope", this->_systAwareGoodJetExpr);
    }

    std::string filterExpr;
    if (this->_systAwareSelectionType == "bJetUpEnvelope")
    {
        if (!rndDS.HasColumn(this->_systAwareBTagBranch))
        {
            rdfWS_utility::messageERROR("SkimControl", "bJetUpEnvelope requested missing btagBranch: " + this->_systAwareBTagBranch);
        }
        rndDS = rndDS.Define(
            "BJetCond_skimJMEUpEnvelope",
            "GoodJetCond_skimJMEUpEnvelope && (" + this->_systAwareBTagBranch + " > " + std::to_string(this->_systAwareBTagThreshold) + ")"
        );
        rndDS = rndDS.Define("nBJet_skimJMEUpEnvelope", "Nonzero(BJetCond_skimJMEUpEnvelope).size()");
        filterExpr = "(nBJet_skimJMEUpEnvelope >= " + std::to_string(this->_systAwareNJet) + ")";
    }
    else
    {
        rndDS = rndDS.Define("nGoodJet_skimJMEUpEnvelope", "Nonzero(GoodJetCond_skimJMEUpEnvelope).size()");
        filterExpr = "(nGoodJet_skimJMEUpEnvelope >= " + std::to_string(this->_systAwareNJet) + ")";
    }

    if (this->_systAwareRequireJVMEnvelope)
        filterExpr += " && JVMweight_skimEnvelope > 0.0";
    rndDS = rndDS.Filter(filterExpr);

    return rndDS;
}

std::vector<std::string> SkimControl::_getBranchArray(ROOT::RDF::RNode rndDS, int isPreliminary)
{
    // keep the branches in the config only and dump into files
    std::vector<std::string> branchArray;
    std::vector<std::string> originalBRs = rndDS.GetColumnNames();
    for (const auto &brName : this->_branchList)
    {
        if (std::find(originalBRs.begin(), originalBRs.end(), brName) == originalBRs.end())
            continue;
        branchArray.push_back(brName);
    }

    return branchArray;
}

////////////////////////////////////////////////// Process the Skim running

void SkimControl::run()
{
    // enhance ctrl+C trigger sensitivity
    instance = this;
    signal(SIGINT, SkimControl::signalHandler);
    signal(SIGTERM, SkimControl::signalHandler);


    auto hasValidEventsTree = [](const std::string &filePath) -> bool
    {
        std::unique_ptr<TFile> fTemp{TFile::Open(filePath.c_str(), "READ")};
        if (!fTemp || fTemp->IsZombie())
        {
            rdfWS_utility::messageWARN("SkimControl", filePath + " not exist!");
            return false;
        }
        if (!fTemp->GetListOfKeys()->FindObject("Events"))
        {
            rdfWS_utility::messageWARN("SkimControl", filePath + " has no Events tree, skip.");
            return false;
        }
        auto eventsTree = dynamic_cast<TTree *>(fTemp->Get("Events"));
        if (!eventsTree)
        {
            rdfWS_utility::messageWARN("SkimControl", filePath + " Events is not a valid TTree, skip.");
            return false;
        }
        return true;
    };

    for (const auto &channel : this->_channels)
    {
        // detect stopping message, handling after this loop
        if (stop_requested)
        {
            rdfWS_utility::messageINFO("SkimControl", "Stop requested. Exiting before processing dataset " + channel);
            break;
        }

        // skip off channels 
        if (!this->_isOn[channel])
            continue;
        
        rdfWS_utility::messageINFO("SkimControl", "Starting processing channel " + channel);
        auto filePaths = this->_samples.value().getFiles(channel);
        if (filePaths.size() == 0)
            continue;
        
        // skim property
        auto isData = this->_isData;
        auto isPreSkim = this->_isPreliminary;

        // in case of MC, for preliminary skimming, need get the MC total genWeight
        double totalGenWeight(0.0);
        // if ((!isData) && isPreSkim)
        //     totalGenWeight = _getTotalGenWeight(filePaths);

        // determine the outDir
        auto outDir = this->_outDir+"/"+ (isData?"data":"mc") +"/"+this->_era+"/"+channel+"/";
        rdfWS_utility::creatingFolder("SkimControl", outDir);

        // enable choose to do merging or not
        if (this->_ifMerging) 
        {
            // collect total genWeight for merged output (MC only, preliminary skim)
            double totalGenWeightMerged = 0.0;
            if (!isData && isPreSkim)
                totalGenWeightMerged = _getTotalGenWeight(filePaths);

            // loading dataframe with chain safety check
            TChain *chDS = new TChain("Events");
            for (const auto &filePath : filePaths)
            {
                std::unique_ptr<TFile> fTemp{ TFile::Open(filePath.c_str(), "READ") };
                if (!fTemp || fTemp->IsZombie())
                {
                    rdfWS_utility::messageWARN("SkimControl", filePath + " not exist!");
                    continue;
                }
                // for the case of non preliminary skimming, I must having the genWeightSum
                if (!isData && !isPreSkim)
                {
                    if (!fTemp->GetListOfKeys()->FindObject("genWeightSum"))
                    {
                        rdfWS_utility::messageERROR("SkimControl", filePath + " has no genWeightSum for non-preliminary skimming. This is very strange, stop running, please check!");
                        return;
                    }
                    auto tempsumhist = dynamic_cast<TH1D*>(fTemp->Get("genWeightSum"));
                    if (!tempsumhist)
                    {
                        rdfWS_utility::messageERROR("SkimControl", filePath + " has not valid genWeightSum for non-preliminary skimming. This is very strange, stop running, please check!");
                        return;
                    }
                    totalGenWeightMerged += tempsumhist->Integral();
                }

                if (!fTemp->GetListOfKeys()->FindObject("Events"))
                {
                    rdfWS_utility::messageWARN("SkimControl", filePath + " has no Events tree, skip chaining.");
                    continue;
                }
                chDS->Add(filePath.c_str());
            }
            if (chDS->GetNtrees() == 0)
            {
                rdfWS_utility::messageWARN("SkimControl", "No valid input files with Events tree remain for channel " + channel + ", skip channel.");
                delete chDS;
                continue;
            }
            ROOT::RDataFrame rdfDS(*chDS);
            ROOT::RDF::RNode rndDS(rdfDS);

            // for preliminary skim, must applying golden json for data
            if (this->_isPreliminary && isData)
                rndDS = this->_preliminaryDeco(rndDS, channel, totalGenWeightMerged);

            // apply the filter
            if (this->_useSystAwareFourJet && !isData)
                rndDS = this->_applySystAwareFourJet(rndDS);
            if (this->_useSystAwareFourJet && !isData && this->_skipNominalFourJetSteps)
            {
                rndDS = this->_skimCut.applyCutSkippingSteps(rndDS, this->_skimSystSkipSteps);
                std::map<std::string, std::string> skimExprReplacements{
                    {"Jet_pt_JEC", "Jet_pt_JEC_skimJMEUpEnvelope"}
                };
                if (this->_systAwareRequireJVMEnvelope)
                    skimExprReplacements["JVMweight"] = "JVMweight_skimEnvelope";
                rndDS = this->_skimCut.applySuffixedCutSubset(
                    rndDS,
                    this->_skimSystSkipSteps,
                    "_skimJMEUpEnvelope",
                    skimExprReplacements);
            }
            else
                rndDS = this->_skimCut.applyCut(rndDS);

            // keep the branches in the config only and dump into files
            auto branchArray = this->_getBranchArray(rndDS, isPreSkim);

            // output
            ROOT::RDF::RSnapshotOptions SSoption;
            SSoption.fCompressionLevel = 6;
            std::string outputPath = outDir + "/" + channel + "_" + this->_skimName + ".root";
            rndDS.Snapshot("Events", outputPath, branchArray, SSoption);

            if (!isData)
            {
                TFile* fRC = new TFile(outputPath.c_str(), "UPDATE");
                auto h_sumw = new TH1D("genWeightSum", "sum of genWeight (this file)", 1, 0.0, 1.0);
                h_sumw->SetBinContent(1, totalGenWeightMerged);
                h_sumw->GetYaxis()->SetTitle("sum(genWeight)");
                h_sumw->Write("");
                fRC->Close();
                delete fRC;
            }

            // collect memory
            delete chDS;
        }
        else
        {
            std::vector<std::pair<std::string, std::string>> jobs;

            for (const auto &filePath: filePaths)
            {
                if (!hasValidEventsTree(filePath))
                {
                    rdfWS_utility::messageWARN("SkimControl", "Skip job construction for invalid input file: " + filePath);
                    continue;
                }
                // file naming structure
                // original: data or mc /era/channel/ (NANOAOD/MINIAOD) or condition/ runNumber / sampleName
                // after Rochester corr, uniformed: data or mc /era/channel/ (runNumber+sample)
                auto outSampleName = filePath.substr(filePath.rfind("/")+1);
                auto dirPart = filePath.substr(0, filePath.rfind("/"));
                auto runNumber =  dirPart.substr(dirPart.rfind("/")+1);
                outSampleName = runNumber+"-"+outSampleName;

                std::string outputPath = outDir + "/" + outSampleName;

                if (std::filesystem::exists(outputPath)) 
                {
                    rdfWS_utility::messageINFO("SkimControl", outputPath + " already exist, skip");
                    continue;
                }
                jobs.emplace_back(filePath, outputPath);
            }

            if (this->_maxFilesPerChannel>0 && jobs.size()>this->_maxFilesPerChannel) 
            {
                jobs.erase(jobs.begin()+this->_maxFilesPerChannel, jobs.end());
            }

            auto totalJobs = jobs.size();
            if (totalJobs==0) continue;
            int batchSize=4; // hard coding now

            // batchSize=1;
            // doing parallel skimming, after each batch, need to write back the weight back
            if (batchSize>1)
            {
                int batches = (totalJobs-1)/batchSize+1;

                // running inside a batch
                for (int i=0; i < batches; i++)
                {
                    rdfWS_utility::messageINFO("SkimControl", "Executing batch "+std::to_string(i));

                    if (stop_requested)
                    {
                        rdfWS_utility::messageINFO("SkimControl", "Stop requested. Exiting before processing dataset " + channel);
                        break;
                    }

                    std::vector<std::unique_ptr<ROOT::RDataFrame>> dfs;

                    using SnapRet_t = ROOT::RDF::RResultPtr<ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager>>;
                    using SumRet_t  = ROOT::RDF::RResultPtr<float>;
                    std::vector<SnapRet_t> snapshots; 
                    std::vector<SumRet_t>  sumws;

                    std::vector<std::string> outPaths; 

                    int nInBatch = batchSize;
                    if (i==batches-1) nInBatch = totalJobs-batchSize*i;

                    dfs.reserve(nInBatch);
                    snapshots.reserve(nInBatch);
                    sumws.reserve(nInBatch);
                    outPaths.reserve(nInBatch);

                    ROOT::RDF::RSnapshotOptions SSoption;
                    SSoption.fCompressionLevel = 4;
                    SSoption.fLazy = true;

                    for (int jobId = batchSize*i; jobId < batchSize*i+nInBatch; jobId++)
                    {
                        auto filePath = jobs[jobId].first;
                        auto outputPath = jobs[jobId].second;

                        // loading dataframe
                        // ROOT::RDataFrame rdfDS("Events", filePath.c_str());
                        rdfWS_utility::messageINFO("SkimControl", "Adding "+filePath);
                        dfs.emplace_back(std::make_unique<ROOT::RDataFrame>("Events", filePath.c_str()));
                        auto& rdfDS = *dfs.back();
                        ROOT::RDF::RNode rndDS(rdfDS);

                        // note could only for data
                        if (!isData)
                            sumws.emplace_back(rndDS.Sum<float>("genWeight"));
                        else
                            sumws.emplace_back(rndDS.Define("dummyWeight", "1.0f").Sum<float>("dummyWeight"));
                            // sumws.emplace_back(rndDS.Count());

                        // for preliminary skim, must applying golden json for data, and rescaling MC weight
                        if (this->_isPreliminary && isData)
                            rndDS = this->_preliminaryDeco(rndDS, channel, totalGenWeight);

                        // apply the filter
                        if (this->_useSystAwareFourJet && !isData)
                            rndDS = this->_applySystAwareFourJet(rndDS);
                        if (this->_useSystAwareFourJet && !isData && this->_skipNominalFourJetSteps)
                        {
                            rndDS = this->_skimCut.applyCutSkippingSteps(rndDS, this->_skimSystSkipSteps);
                            std::map<std::string, std::string> skimExprReplacements{
                                {"Jet_pt_JEC", "Jet_pt_JEC_skimJMEUpEnvelope"}
                            };
                            if (this->_systAwareRequireJVMEnvelope)
                                skimExprReplacements["JVMweight"] = "JVMweight_skimEnvelope";
                            rndDS = this->_skimCut.applySuffixedCutSubset(
                                rndDS,
                                this->_skimSystSkipSteps,
                                "_skimJMEUpEnvelope",
                                skimExprReplacements);
                        }
                        else
                            rndDS = this->_skimCut.applyCut(rndDS);

                        // keep the branches in the config only and dump into files
                        auto branchArray = this->_getBranchArray(rndDS, isPreSkim);

                        // output
                        auto snap = rndDS.Snapshot("Events", outputPath, branchArray, SSoption);
                        outPaths.emplace_back(outputPath);
                        snapshots.push_back(snap);
                    }
                    std::vector<ROOT::RDF::RResultHandle> handles;
                    handles.reserve(snapshots.size() + sumws.size());
                    for (auto &s : snapshots) handles.emplace_back(ROOT::RDF::RResultHandle(s));
                    for (auto &w : sumws)     handles.emplace_back(ROOT::RDF::RResultHandle(w));
                    unsigned ran = ROOT::RDF::RunGraphs(handles);

                    // after a batch running, padding the total weight
                    for (int idx_inB = 0; idx_inB < nInBatch; idx_inB++)
                    {
                        auto outputPath = outPaths[idx_inB];
                        const double totalWeight = sumws[idx_inB].GetValue();

                        TFile* fRC = new TFile(outputPath.c_str(), "UPDATE");
                        auto h_sumw = new TH1D("genWeightSum", "sum of genWeight (this file)", 1, 0.0, 1.0);
                        h_sumw->SetBinContent(1, totalWeight);
                        h_sumw->GetYaxis()->SetTitle("sum(genWeight)");
                        h_sumw->Write("");
                        fRC->Close();
                        delete fRC;
                    }
                }
            }
            else
            {
                rdfWS_utility::messageINFO("SkimControl", "Executing non batch ");

                std::vector<std::unique_ptr<ROOT::RDataFrame>> dfs;

                using SnapRet_t = ROOT::RDF::RResultPtr<ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager>>;
                std::vector<SnapRet_t> snapshots; 

                ROOT::RDF::RSnapshotOptions SSoption;
                SSoption.fCompressionLevel = 4;
                SSoption.fLazy = false;

                for (int jobId = 0; jobId < totalJobs; jobId++)
                {
                    auto filePath = jobs[jobId].first;
                    auto outputPath = jobs[jobId].second;

                    if (stop_requested)
                    {
                        rdfWS_utility::messageINFO("SkimControl", "Stop requested. Exiting before processing dataset " + channel);
                        break;
                    }

                    // loading dataframe
                    // ROOT::RDataFrame rdfDS("Events", filePath.c_str());
                    rdfWS_utility::messageINFO("SkimControl", "Adding "+filePath);
                    dfs.emplace_back(std::make_unique<ROOT::RDataFrame>("Events", filePath.c_str()));
                    auto& rdfDS = *dfs.back();
                    ROOT::RDF::RNode rndDS(rdfDS);

                    // for preliminary skim, must applying golden json for data, and rescaling MC weight
                    if (this->_isPreliminary)
                        rndDS = this->_preliminaryDeco(rndDS, channel, totalGenWeight);

                    // apply the filter
                    if (this->_useSystAwareFourJet && !isData)
                        rndDS = this->_applySystAwareFourJet(rndDS);
                    if (this->_useSystAwareFourJet && !isData && this->_skipNominalFourJetSteps)
                    {
                        rndDS = this->_skimCut.applyCutSkippingSteps(rndDS, this->_skimSystSkipSteps);
                        std::map<std::string, std::string> skimExprReplacements{
                            {"Jet_pt_JEC", "Jet_pt_JEC_skimJMEUpEnvelope"}
                        };
                        if (this->_systAwareRequireJVMEnvelope)
                            skimExprReplacements["JVMweight"] = "JVMweight_skimEnvelope";
                        rndDS = this->_skimCut.applySuffixedCutSubset(
                            rndDS,
                            this->_skimSystSkipSteps,
                            "_skimJMEUpEnvelope",
                            skimExprReplacements);
                    }
                    else
                        rndDS = this->_skimCut.applyCut(rndDS);

                    // keep the branches in the config only and dump into files
                    auto branchArray = this->_getBranchArray(rndDS, isPreSkim);

                    // output
                    rndDS.Snapshot("Events", outputPath, branchArray, SSoption);
                }
            }
        }
    }
}
