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
    ROOT::VecOps::RVec<Float_t> Jet_neHEF;
    ROOT::VecOps::RVec<Float_t> Jet_neEmEF;
    ROOT::VecOps::RVec<UChar_t> Jet_nConstituents;
    ROOT::VecOps::RVec<Float_t> Jet_muEF;
    ROOT::VecOps::RVec<Float_t> Jet_chHEF;
    ROOT::VecOps::RVec<UChar_t> Jet_chMultiplicity;
    ROOT::VecOps::RVec<Float_t> Jet_chEmEF;
    ROOT::VecOps::RVec<UChar_t> Jet_neMultiplicity;
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

    auto v_neHEF = dfs.Take<ROOT::VecOps::RVec<float>>("Jet_neHEF");
    auto v_neEm = dfs.Take<ROOT::VecOps::RVec<float>>("Jet_neEmEF");
    auto v_nC = dfs.Take<ROOT::VecOps::RVec<UChar_t>>("Jet_nConstituents");
    auto v_muEF = dfs.Take<ROOT::VecOps::RVec<float>>("Jet_muEF");
    auto v_chHEF = dfs.Take<ROOT::VecOps::RVec<float>>("Jet_chHEF");
    auto v_chM = dfs.Take<ROOT::VecOps::RVec<UChar_t>>("Jet_chMultiplicity");
    auto v_chEm = dfs.Take<ROOT::VecOps::RVec<float>>("Jet_chEmEF");
    auto v_neM = dfs.Take<ROOT::VecOps::RVec<UChar_t>>("Jet_neMultiplicity");

    const auto N = runs->size();
    std::unordered_map<EvtKey, Packed, EvtKeyHash> table;
    table.reserve(N);

    for (std::size_t i = 0; i < N; ++i)
    {
        EvtKey key{(*runs)[i], (*lumis)[i], (*events)[i]};
        Packed val{
            (*v_neHEF)[i],
            (*v_neEm)[i],
            (*v_nC)[i],
            (*v_muEF)[i],
            (*v_chHEF)[i],
            (*v_chM)[i],
            (*v_chEm)[i],
            (*v_neM)[i]};
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

    auto df2 = dft.Define("Jet_neHEF", [fetchF](unsigned int r, unsigned int l, ULong64_t e)
                          {
            if (auto p = fetchF(r,l,e)) return p->Jet_neHEF;
            return ROOT::VecOps::RVec<float>{}; }, {"run", "luminosityBlock", "event"})
                   .Define("Jet_neEmEF", [fetchF](unsigned int r, unsigned int l, ULong64_t e)
                           {
            if (auto p = fetchF(r,l,e)) return p->Jet_neEmEF;
            return ROOT::VecOps::RVec<float>{}; }, {"run", "luminosityBlock", "event"})
                   .Define("Jet_nConstituents", [fetchF](unsigned int r, unsigned int l, ULong64_t e)
                           {
            if (auto p = fetchF(r,l,e)) return p->Jet_nConstituents;
            return ROOT::VecOps::RVec<unsigned char>{}; }, {"run", "luminosityBlock", "event"})
                   .Define("Jet_muEF", [fetchF](unsigned int r, unsigned int l, ULong64_t e)
                           {
            if (auto p = fetchF(r,l,e)) return p->Jet_muEF;
            return ROOT::VecOps::RVec<float>{}; }, {"run", "luminosityBlock", "event"})
                   .Define("Jet_chHEF", [fetchF](unsigned int r, unsigned int l, ULong64_t e)
                           {
            if (auto p = fetchF(r,l,e)) return p->Jet_chHEF;
            return ROOT::VecOps::RVec<float>{}; }, {"run", "luminosityBlock", "event"})
                   .Define("Jet_chMultiplicity", [fetchF](unsigned int r, unsigned int l, ULong64_t e)
                           {
            if (auto p = fetchF(r,l,e)) return p->Jet_chMultiplicity;
            return ROOT::VecOps::RVec<unsigned char>{}; }, {"run", "luminosityBlock", "event"})
                   .Define("Jet_chEmEF", [fetchF](unsigned int r, unsigned int l, ULong64_t e)
                           {
            if (auto p = fetchF(r,l,e)) return p->Jet_chEmEF;
            return ROOT::VecOps::RVec<float>{}; }, {"run", "luminosityBlock", "event"})
                   .Define("Jet_neMultiplicity", [fetchF](unsigned int r, unsigned int l, ULong64_t e)
                           {
            if (auto p = fetchF(r,l,e)) return p->Jet_neMultiplicity;
            return ROOT::VecOps::RVec<unsigned char>{}; }, {"run", "luminosityBlock", "event"});

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
    SampleControl samples("json/samples/RunIII2024Summer24NanoAODv15_MuonCorr.json");
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

        std::unordered_map<std::string, std::string> channelTable{
            {"DY2L", "/data2/common/NanoAOD/mc/v15/RunIII2024Summer24NanoAODv15/DYto2Mu-2Jets_Bin-MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v6/"},
            {"DY2L_low", "/data2/common/NanoAOD/mc/v15/RunIII2024Summer24NanoAODv15/DYto2Mu_Bin-MLL-10to50_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/"},
            {"ttbarDL", "/data2/common/NanoAOD/mc/v15/RunIII2024Summer24NanoAODv15/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v3/"},
            {"ttbarSL", "/data2/common/NanoAOD/mc/v15/RunIII2024Summer24NanoAODv15/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/"},
            {"ttbar4Q", "/data2/common/NanoAOD/mc/v15/RunIII2024Summer24NanoAODv15/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/"},
            {"WW2L2Nu", "/data2/common/NanoAOD/mc/v15/RunIII2024Summer24NanoAODv15/WWto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/"},
            {"WZ", "/data2/common/NanoAOD/mc/v15/RunIII2024Summer24NanoAODv15/WZ_TuneCP5_13p6TeV_pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/"},
            {"ZZ", "/data2/common/NanoAOD/mc/v15/RunIII2024Summer24NanoAODv15/ZZ_TuneCP5_13p6TeV_pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/"},
            {"WWW", "/data2/common/NanoAOD/mc/v15/RunIII2024Summer24NanoAODv15/WWW-4F_TuneCP5_13p6TeV_amcatnlo-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/"},
            {"WWZ", "/data2/common/NanoAOD/mc/v15/RunIII2024Summer24NanoAODv15/WWZ-4F_TuneCP5_13p6TeV_amcatnlo-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/"},
            {"WZZ", "/data2/common/NanoAOD/mc/v15/RunIII2024Summer24NanoAODv15/WZZ-5F_TuneCP5_13p6TeV_amcatnlo-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/"},
            {"ZZZ", "/data2/common/NanoAOD/mc/v15/RunIII2024Summer24NanoAODv15/ZZZ-5F_TuneCP5_13p6TeV_amcatnlo-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/"},
            {"TWm2L", "/data2/common/NanoAOD/mc/v15/RunIII2024Summer24NanoAODv15/TWminusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/"},
            {"TbarWp2L", "/data2/common/NanoAOD/mc/v15/RunIII2024Summer24NanoAODv15/TbarWplusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/"},
            {"WZ2L2Q", "/data2/common/NanoAOD/mc/v15/RunIII2024Summer24NanoAODv15/WZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/"},
            {"ZZ2L2Q", "/data2/common/NanoAOD/mc/v15/RunIII2024Summer24NanoAODv15/ZZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/"},
            {"ZZ2L2Nu", "/data2/common/NanoAOD/mc/v15/RunIII2024Summer24NanoAODv15/ZZto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/"},
            {"WZ3l", "/data2/common/NanoAOD/mc/v15/RunIII2024Summer24NanoAODv15/WZto3LNu_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/"},
            {"ZZ4l", "/data2/common/NanoAOD/mc/v15/RunIII2024Summer24NanoAODv15/ZZto4L_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/"},
            {"TTWW", "/data2/common/NanoAOD/mc/v15/RunIII2024Summer24NanoAODv15/TTWW_TuneCP5_13p6TeV_madgraph-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/"},
            {"TTWZ", "/data2/common/NanoAOD/mc/v15/RunIII2024Summer24NanoAODv15/TTWZ_TuneCP5_13p6TeV_madgraph-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/"},
            {"TTHBB", "/data2/common/NanoAOD/mc/v15/RunIII2024Summer24NanoAODv15/TTH-Hto2B_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/"},
            {"TTHnonBB", "/data2/common/NanoAOD/mc/v15/RunIII2024Summer24NanoAODv15/TTH-HtoNon2B_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/"},
            {"QCD_80_120_mu", "/data2/common/NanoAOD/mc/v15/RunIII2024Summer24NanoAODv15/QCD_Bin-PT-80to120_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/"},
            {"QCD_120_170_mu", "/data2/common/NanoAOD/mc/v15/RunIII2024Summer24NanoAODv15/QCD_Bin-PT-120to170_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/"},
            {"QCD_170_300_mu", "/data2/common/NanoAOD/mc/v15/RunIII2024Summer24NanoAODv15/QCD_Bin-PT-170to300_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/"},
        };

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

        run(srcFile, tgtFile, outFile);
    }
    return 0;
}
