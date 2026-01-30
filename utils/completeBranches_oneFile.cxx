// make_muon_friend.C
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TFile.h>
#include <TTree.h>
#include <unordered_set>
#include <string>
#include <map>
#include <filesystem>
#include <algorithm>
#include <iostream>
#include <TH1.h>

#include "SampleControl.h"
#include "Utility.h"


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


bool process(const std::string &skimFile,
             const std::string &rawFile,
             const std::string &branchJson,
             const std::string &outFriend)
{
    const std::filesystem::path outPath(outFriend);
    if (!outPath.parent_path().empty() && !std::filesystem::exists(outPath.parent_path()))
    {
        std::filesystem::create_directories(outPath.parent_path());
    }
    if (std::filesystem::exists(outFriend))
    {
        TFile fout(outFriend.c_str(), "READ");
        if (fout.IsZombie())
        {
            rdfWS_utility::messageWARN("completeBranches_oneFile", "Output is zombie, removing: " + outFriend);
            fout.Close();
            std::filesystem::remove(outFriend);
        }
        else
        {
            rdfWS_utility::messageWARN("completeBranches_oneFile", "Output exists, skip: " + outFriend);
            return true;
        }
    }

    rdfWS_utility::JsonObject branchObj(rdfWS_utility::readJson("completeBranches", branchJson), "Branch Config");
    auto branches = branchObj.get<std::vector<std::string>>();
    if (std::find(branches.begin(), branches.end(), "genWeightSum") == branches.end())
    {
        branches.push_back("genWeightSum");
    }

    // 1) build key set from skim (only run/lumi/event)
    ROOT::RDataFrame dskim("Events", skimFile);
    auto runs = dskim.Take<unsigned int>("run");
    auto lumis = dskim.Take<unsigned int>("luminosityBlock");
    auto evts = dskim.Take<ULong64_t>("event");

    std::unordered_set<EvtKey, EvtKeyHash> keys;
    keys.reserve(runs->size());
    for (size_t i = 0; i < runs->size(); ++i)
    {
        keys.insert(EvtKey{(*runs)[i], (*lumis)[i], (*evts)[i]});
    }

    // 2) scan raw, keep only events present in skim
    ROOT::RDataFrame draw("Events", rawFile);
    auto *pkeys = &keys;

    std::vector<std::string> savedBr;
    std::vector<std::string> fileBr = draw.GetColumnNames();
    for (const auto &brName : branches)
    {
        if (std::find(fileBr.begin(), fileBr.end(), brName) == fileBr.end())
            continue;
        savedBr.push_back(brName);
    }

    auto df = draw.Filter([pkeys](unsigned int r, unsigned int l, ULong64_t e){ return pkeys->find(EvtKey{r, l, e}) != pkeys->end(); }, {"run", "luminosityBlock", "event"});
    df.Snapshot("Events", outFriend, savedBr);

    // copy genWeightSum histogram from skim file (outside Events tree)
    {
        TFile fskim(skimFile.c_str(), "READ");
        if (!fskim.IsOpen())
        {
            rdfWS_utility::messageWARN("completeBranches_oneFile", "Failed to open skim file for genWeightSum: " + skimFile);
            return false;
        }
        auto *h = dynamic_cast<TH1*>(fskim.Get("genWeightSum"));
        if (!h)
        {
            rdfWS_utility::messageERROR("completeBranches_oneFile", "genWeightSum not found in skim file: " + skimFile);
            std::filesystem::remove(outFriend);
            return false;
        }
        TFile fout(outFriend.c_str(), "UPDATE");
        if (!fout.IsOpen())
        {
            rdfWS_utility::messageWARN("completeBranches_oneFile", "Failed to open output file for genWeightSum: " + outFriend);
            return false;
        }
        h->SetDirectory(&fout);
        h->Write("genWeightSum", TObject::kOverwrite);
    }
    return true;
}

int main(int argc, char *argv[])
{
    std::string skimmedFile, rawFile, branchConfig, outFile;
    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        if (arg == "--skimmed")
        {
            if (i + 1 >= argc)
            {
                rdfWS_utility::messageERROR("completeBranches_oneFile", "Missing value for --skimmed");
                return 1;
            }
            skimmedFile = argv[++i];
        }
        else if (arg == "--rawFile")
        {
            if (i + 1 >= argc)
            {
                rdfWS_utility::messageERROR("completeBranches_oneFile", "Missing value for --rawFile");
                return 1;
            }
            rawFile = argv[++i];
        }
        else if (arg == "--branches")
        {
            if (i + 1 >= argc)
            {
                rdfWS_utility::messageERROR("completeBranches_oneFile", "Missing value for --branches");
                return 1;
            }
            branchConfig = argv[++i];
        }
        else if (arg == "--outFile")
        {
            if (i + 1 >= argc)
            {
                rdfWS_utility::messageERROR("completeBranches_oneFile", "Missing value for --outFile");
                return 1;
            }
            outFile = argv[++i];
        }
        else
        {
            rdfWS_utility::messageWARN("completeBranches_oneFile",
                                       "Unknown argument: " + arg + ". Expected: --skimmed --rawFile --branches --outFile");
            return 1;
        }
    }

    if (skimmedFile.empty() || rawFile.empty() || branchConfig.empty() || outFile.empty())
    {
        rdfWS_utility::messageERROR(
            "completeBranches_oneFile",
            "Input template wrong. Example: \"completeBranches_oneFile --skimmed <skimmed file> --rawFile <raw file> --branches <branch json> --outFile <output file>\"");
        return 1;
    }

    if (!process(skimmedFile, rawFile, branchConfig, outFile))
        return 1;
    return 0;
}
