// main.cpp
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TROOT.h>
#include <filesystem>
#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>
#include <tuple>
#include <chrono>
#include <iomanip>
#include <ctime>

#include "SampleControl.h"

// file-wise memory storage
//////////////////////////////////////////////////
// event identifier
struct EvtKey
{
    unsigned int run;
    unsigned int lumi;
    ULong64_t event;
    bool operator==(const EvtKey &o) const
    {
        return run == o.run && lumi == o.lumi && event == o.event;
    }
};
struct EvtKeyHash
{
    std::size_t operator()(const EvtKey &k) const noexcept
    {
        auto h1 = std::hash<unsigned int>{}(k.run);
        auto h2 = std::hash<unsigned int>{}(k.lumi);
        auto h3 = std::hash<ULong64_t>{}(k.event);
        return ((h1 * 1315423911u) ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2))) ^ h3;
    }
};

// hash map target
struct Packed
{
    ROOT::VecOps::RVec<Float_t> Jet_rawFactor;
    Float_t Rho_fixedGridRhoFastjetAll;
    ROOT::VecOps::RVec<Short_t> Jet_genJetIdx;
    Int_t nGenJet;
    ROOT::VecOps::RVec<Float_t> GenJet_pt;
    ROOT::VecOps::RVec<Float_t> GenJet_eta;
    ROOT::VecOps::RVec<Float_t> GenJet_phi;
    Float_t RawPuppiMET_phi;
    Float_t RawPuppiMET_pt;
    Float_t RawPuppiMET_sumEt;
    ROOT::VecOps::RVec<Float_t> Jet_muonSubtrFactor;
};

struct DataPacked
{
    ROOT::VecOps::RVec<Float_t> Jet_rawFactor;
    Float_t Rho_fixedGridRhoFastjetAll;
    Float_t RawPuppiMET_phi;
    Float_t RawPuppiMET_pt;
    Float_t RawPuppiMET_sumEt;
    ROOT::VecOps::RVec<Float_t> Jet_muonSubtrFactor;
};

// operation function
//////////////////////////////////////////////////
static std::string now_str()
{
    using clock = std::chrono::system_clock;
    auto now = clock::now();
    std::time_t t = clock::to_time_t(now);
    std::tm tm{};
    localtime_r(&t, &tm);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S") << " ";
    return oss.str();
}

void run(std::string srcFile, std::string tgtFile, std::string outFile)
{
    // read into memory
    std::cout << now_str() << "staring loading from the original: " << srcFile << std::endl;
    ROOT::RDataFrame dfs("Events", srcFile);

    auto runs = dfs.Take<unsigned int>("run");
    auto lumis = dfs.Take<unsigned int>("luminosityBlock");
    auto events = dfs.Take<ULong64_t>("event");

    auto v_Jet_rawFactor = dfs.Take<ROOT::VecOps::RVec<Float_t>>("Jet_rawFactor");
    auto v_Rho_fixedGridRhoFastjetAll = dfs.Take<Float_t>("Rho_fixedGridRhoFastjetAll");
    auto v_Jet_genJetIdx = dfs.Take<ROOT::VecOps::RVec<Short_t>>("Jet_genJetIdx");
    auto v_nGenJet = dfs.Take<Int_t>("nGenJet");
    auto v_GenJet_pt = dfs.Take<ROOT::VecOps::RVec<Float_t>>("GenJet_pt");
    auto v_GenJet_eta = dfs.Take<ROOT::VecOps::RVec<Float_t>>("GenJet_eta");
    auto v_GenJet_phi = dfs.Take<ROOT::VecOps::RVec<Float_t>>("GenJet_phi");
    auto v_RawPuppiMET_phi = dfs.Take<Float_t>("RawPuppiMET_phi");
    auto v_RawPuppiMET_pt = dfs.Take<Float_t>("RawPuppiMET_pt");
    auto v_RawPuppiMET_sumEt = dfs.Take<Float_t>("RawPuppiMET_sumEt");
    auto v_Jet_muonSubtrFactor = dfs.Take<ROOT::VecOps::RVec<Float_t>>("Jet_muonSubtrFactor");

    const auto N = runs->size();
    std::unordered_map<EvtKey, Packed, EvtKeyHash> table;
    table.reserve(N);

    for (std::size_t i = 0; i < N; ++i)
    {
        EvtKey key{(*runs)[i], (*lumis)[i], (*events)[i]};
        Packed val{
            (*v_Jet_rawFactor)[1],
            (*v_Rho_fixedGridRhoFastjetAll)[1],
            (*v_Jet_genJetIdx)[1],
            (*v_nGenJet)[1],
            (*v_GenJet_pt)[1],
            (*v_GenJet_eta)[1],
            (*v_GenJet_phi)[1],
            (*v_RawPuppiMET_phi)[1],
            (*v_RawPuppiMET_pt)[1],
            (*v_RawPuppiMET_sumEt)[1],
            (*v_Jet_muonSubtrFactor)[1]};
        table.emplace(std::move(key), std::move(val));
    }
    std::cout << now_str() << "Loaded " << table.size() << " events into memory.\n";

    // output into the target file
    ROOT::RDataFrame dft("Events", tgtFile);

    auto *pTable = &table;
    auto fetchF = [pTable](unsigned int r, unsigned int l, ULong64_t e) -> const Packed *
    {
        auto it = pTable->find(EvtKey{r, l, e});
        return (it == pTable->end()) ? nullptr : &(it->second);
    };


    auto df2 =  dft.Define("Jet_rawFactor", [fetchF](unsigned int r, unsigned int l, ULong64_t e) {
            if (auto p = fetchF(r,l,e)) return p->Jet_rawFactor;
            return ROOT::VecOps::RVec<Float_t>{}; }, {"run", "luminosityBlock", "event"})
                   .Define("Rho_fixedGridRhoFastjetAll", [fetchF](unsigned int r, unsigned int l, ULong64_t e) {
            if (auto p = fetchF(r,l,e)) return p->Rho_fixedGridRhoFastjetAll;
            return Float_t(-999.); }, {"run", "luminosityBlock", "event"})
                   .Define("Jet_genJetIdx", [fetchF](unsigned int r, unsigned int l, ULong64_t e) {
            if (auto p = fetchF(r,l,e)) return p->Jet_genJetIdx;
            return ROOT::VecOps::RVec<Short_t>{}; }, {"run", "luminosityBlock", "event"})
                   .Define("nGenJet", [fetchF](unsigned int r, unsigned int l, ULong64_t e) {
            if (auto p = fetchF(r,l,e)) return p->nGenJet;
            return 0; }, {"run", "luminosityBlock", "event"})
                   .Define("GenJet_pt", [fetchF](unsigned int r, unsigned int l, ULong64_t e) {
            if (auto p = fetchF(r,l,e)) return p->GenJet_pt;
            return ROOT::VecOps::RVec<Float_t>{}; }, {"run", "luminosityBlock", "event"})
                   .Define("GenJet_eta", [fetchF](unsigned int r, unsigned int l, ULong64_t e) {
            if (auto p = fetchF(r,l,e)) return p->GenJet_eta;
            return ROOT::VecOps::RVec<Float_t>{}; }, {"run", "luminosityBlock", "event"})
                   .Define("GenJet_phi", [fetchF](unsigned int r, unsigned int l, ULong64_t e) {
            if (auto p = fetchF(r,l,e)) return p->GenJet_phi;
            return ROOT::VecOps::RVec<Float_t>{}; }, {"run", "luminosityBlock", "event"})
                   .Define("RawPuppiMET_phi", [fetchF](unsigned int r, unsigned int l, ULong64_t e) {
            if (auto p = fetchF(r,l,e)) return p->RawPuppiMET_phi;
            return Float_t(-999.); }, {"run", "luminosityBlock", "event"})
                   .Define("RawPuppiMET_pt", [fetchF](unsigned int r, unsigned int l, ULong64_t e) {
            if (auto p = fetchF(r,l,e)) return p->RawPuppiMET_pt;
            return Float_t(-999.); }, {"run", "luminosityBlock", "event"})
                   .Define("RawPuppiMET_sumEt", [fetchF](unsigned int r, unsigned int l, ULong64_t e) {
            if (auto p = fetchF(r,l,e)) return p->RawPuppiMET_sumEt;
            return Float_t(-999.); }, {"run", "luminosityBlock", "event"})
                   .Define("Jet_muonSubtrFactor", [fetchF](unsigned int r, unsigned int l, ULong64_t e) {
            if (auto p = fetchF(r,l,e)) return p->Jet_muonSubtrFactor;
            return ROOT::VecOps::RVec<Float_t>{}; }, {"run", "luminosityBlock", "event"});

    // output with column increment
    df2.Snapshot("Events", outFile);

    std::cout << now_str() << "Wrote: " << outFile << "\n";
}


void runData(std::string srcFile, std::string tgtFile, std::string outFile)
{
    // read into memory
    std::cout << now_str() << "staring loading from the original: " << srcFile << std::endl;
    ROOT::RDataFrame dfs("Events", srcFile);

    auto runs = dfs.Take<unsigned int>("run");
    auto lumis = dfs.Take<unsigned int>("luminosityBlock");
    auto events = dfs.Take<ULong64_t>("event");

    auto v_Jet_rawFactor = dfs.Take<ROOT::VecOps::RVec<Float_t>>("Jet_rawFactor");
    auto v_Rho_fixedGridRhoFastjetAll = dfs.Take<Float_t>("Rho_fixedGridRhoFastjetAll");
    auto v_RawPuppiMET_phi = dfs.Take<Float_t>("RawPuppiMET_phi");
    auto v_RawPuppiMET_pt = dfs.Take<Float_t>("RawPuppiMET_pt");
    auto v_RawPuppiMET_sumEt = dfs.Take<Float_t>("RawPuppiMET_sumEt");
    auto v_Jet_muonSubtrFactor = dfs.Take<ROOT::VecOps::RVec<Float_t>>("Jet_muonSubtrFactor");

    const auto N = runs->size();
    std::unordered_map<EvtKey, DataPacked, EvtKeyHash> table;
    table.reserve(N);

    for (std::size_t i = 0; i < N; ++i)
    {
        EvtKey key{(*runs)[i], (*lumis)[i], (*events)[i]};
        DataPacked val{
            (*v_Jet_rawFactor)[1],
            (*v_Rho_fixedGridRhoFastjetAll)[1],
            (*v_RawPuppiMET_phi)[1],
            (*v_RawPuppiMET_pt)[1],
            (*v_RawPuppiMET_sumEt)[1],
            (*v_Jet_muonSubtrFactor)[1]};
        table.emplace(std::move(key), std::move(val));
    }
    std::cout << now_str() << "Loaded " << table.size() << " events into memory.\n";

    // output into the target file
    ROOT::RDataFrame dft("Events", tgtFile);

    auto *pTable = &table;
    auto fetchF = [pTable](unsigned int r, unsigned int l, ULong64_t e) -> const DataPacked *
    {
        auto it = pTable->find(EvtKey{r, l, e});
        return (it == pTable->end()) ? nullptr : &(it->second);
    };


    auto df2 =  dft.Define("Jet_rawFactor", [fetchF](unsigned int r, unsigned int l, ULong64_t e) {
            if (auto p = fetchF(r,l,e)) return p->Jet_rawFactor;
            return ROOT::VecOps::RVec<Float_t>{}; }, {"run", "luminosityBlock", "event"})
                   .Define("Rho_fixedGridRhoFastjetAll", [fetchF](unsigned int r, unsigned int l, ULong64_t e) {
            if (auto p = fetchF(r,l,e)) return p->Rho_fixedGridRhoFastjetAll;
            return Float_t(-999.); }, {"run", "luminosityBlock", "event"})
                   .Define("RawPuppiMET_phi", [fetchF](unsigned int r, unsigned int l, ULong64_t e) {
            if (auto p = fetchF(r,l,e)) return p->RawPuppiMET_phi;
            return Float_t(-999.); }, {"run", "luminosityBlock", "event"})
                   .Define("RawPuppiMET_pt", [fetchF](unsigned int r, unsigned int l, ULong64_t e) {
            if (auto p = fetchF(r,l,e)) return p->RawPuppiMET_pt;
            return Float_t(-999.); }, {"run", "luminosityBlock", "event"})
                   .Define("RawPuppiMET_sumEt", [fetchF](unsigned int r, unsigned int l, ULong64_t e) {
            if (auto p = fetchF(r,l,e)) return p->RawPuppiMET_sumEt;
            return Float_t(-999.); }, {"run", "luminosityBlock", "event"})
                   .Define("Jet_muonSubtrFactor", [fetchF](unsigned int r, unsigned int l, ULong64_t e) {
            if (auto p = fetchF(r,l,e)) return p->Jet_muonSubtrFactor;
            return ROOT::VecOps::RVec<Float_t>{}; }, {"run", "luminosityBlock", "event"});

    // output with column increment
    df2.Snapshot("Events", outFile);

    std::cout << now_str() << "Wrote: " << outFile << "\n";
}


std::vector<std::string> getDataFiles(std::string era, std::string channel)
{
    SampleControl samples("json/samples/" + era + "_forJVM.json");
    return samples.getFiles(channel);
}

std::vector<std::string> getMCFiles(std::string channel)
{
    SampleControl samples("json/samples/RunIII2024Summer24NanoAODv15_forJVM.json");
    return samples.getFiles(channel);
}

std::string getSource(const std::string &targetFile, const std::string &era, const std::string &channel, bool isData)
{
    if (isData)
    {
        // e.g
        // /data2/common/skimmed_NanoAOD/skim_1027_MVAtight_TTHH/data/Run2024C/Muon0/2530000-a1c21041-2199-4745-bb97-1da059203510_Rcorr.root
        // to /data2/common/NanoAOD/data/v15/Run2024C/Muon0/NANOAOD/MINIv6NANOv15-v1/2530000/a1c21041-2199-4745-bb97-1da059203510.root
        auto tFileName = targetFile.substr(targetFile.rfind("/") + 1);
        auto runNum = tFileName.substr(0, tFileName.find("-"));
        auto oFileName = tFileName.substr(tFileName.find("-") + 1);
        oFileName = oFileName.substr(0, oFileName.find("_Rcorr")) + ".root";
        // special cases:
        // Run2024G and Run2024H have Muon1 with MINIv6NANOv15-v1
        // Run2024I has difference in the Muon0 and Muon0_add
        if ((era=="Run2024G" || era == "Run2024H") && channel =="Muon1")
        {
            return "/data2/common/NanoAOD/data/v15/" + era + "/" + channel + "/NANOAOD/MINIv6NANOv15-v2/" + runNum + "/" + oFileName;
        }
        if (era=="Run2024I")
        {
            // checking the last forlder
            auto tFileDir = targetFile.substr(0, targetFile.rfind("/"));
            auto subChannel = tFileDir.substr(tFileDir.rfind("/")+1);
            if (subChannel.find("_") != std::string::npos)
            {
                return "/data2/common/NanoAOD/data/v15/" + era + "/" + channel + "/NANOAOD/MINIv6NANOv15_v2-v1/" + runNum + "/" + oFileName;
            }
            return "/data2/common/NanoAOD/data/v15/" + era + "/" + channel + "/NANOAOD/MINIv6NANOv15-v1/" + runNum + "/" + oFileName;
        }
        // normal case
        return "/data2/common/NanoAOD/data/v15/" + era + "/" + channel + "/NANOAOD/MINIv6NANOv15-v1/" + runNum + "/" + oFileName;
    }
    else
    {
        auto tFileName = targetFile.substr(targetFile.rfind("/") + 1);
        auto runNum = tFileName.substr(0, tFileName.find("-"));
        auto oFileName = tFileName.substr(tFileName.find("-") + 1);
        oFileName = oFileName.substr(0, oFileName.find("_Rcorr")) + ".root";

        SampleControl samples("json/samples/" + era + ".json");
        auto channelTable = samples.getAllDirs();

        return channelTable.at(channel) + runNum + "/" + oFileName;
    }
}

std::string getOutput(const std::string& targetFile)
{
    // e.g. /data2/common/skimmed_NanoAOD/skim_1027_MVAtight_TTHH/mc/RunIII2024Summer24NanoAODv15/DY2L/2520000-a1139b21-2ecd-4e40-8fc0-1fddccf89830_Rcorr.root
    // TO /data2/common/skimmed_NanoAOD/skim_1027_MVAtight_TTHH/mc/RunIII2024Summer24NanoAODv15/DY2L_corr/2520000-a1139b21-2ecd-4e40-8fc0-1fddccf89830_Rcorr.root
    auto tFileName = targetFile.substr(targetFile.rfind("/") + 1);
    auto tFileDir = targetFile.substr(0, targetFile.rfind("/"));
    if (!std::filesystem::exists(tFileDir+"_corr/"))
    {
        std::filesystem::create_directory(tFileDir+"_corr/");
    }
    return tFileDir+"_corr/"+tFileName;
}

// main function
//////////////////////////////////////////////////
int main(int argc, char **argv)
{
    if (argc < 3)
    {
        std::cout << "Need providing period and channels" << std::endl;
        std::cout << "Note only for 2024 now as a quick fix" << std::endl;
        exit(1);
    }
    std::string era = argv[1];
    std::string channel = argv[2];
    bool isData = (era != "RunIII2024Summer24NanoAODv15");

    ROOT::EnableImplicitMT();

    // get files for processing
    std::vector<std::string> allTargetFiles;
    if (isData)
        allTargetFiles = getDataFiles(era, channel);
    else
        allTargetFiles = getMCFiles(channel);

    for (const auto& tgtFile: allTargetFiles)
    {
        auto srcFile = getSource(tgtFile, era, channel, isData);
        auto outFile = getOutput(tgtFile);

        if (std::filesystem::exists(outFile))
            continue;

        if (isData) 
            runData(srcFile, tgtFile, outFile);
        else
            run(srcFile, tgtFile, outFile);
    }
    return 0;
}
