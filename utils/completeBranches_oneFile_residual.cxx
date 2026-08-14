// make_muon_friend.C
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <string>
#include <filesystem>
#include <vector>

#include "SampleControl.h"
#include "Utility.h"


bool process(const std::string &skimFile,
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

    // preload genWeightSum from skim file (required)
    TFile fskim(skimFile.c_str(), "READ");
    if (!fskim.IsOpen())
    {
        rdfWS_utility::messageWARN("completeBranches_oneFile", "Failed to open skim file: " + skimFile);
        return false;
    }
    auto *h = dynamic_cast<TH1 *>(fskim.Get("genWeightSum"));
    if (!h)
    {
        rdfWS_utility::messageERROR("completeBranches_oneFile", "genWeightSum not found in skim file: " + skimFile);
        return false;
    }
    TH1 *hClone = dynamic_cast<TH1 *>(h->Clone("genWeightSum"));
    hClone->SetDirectory(nullptr);
    TTree *skimEvents = dynamic_cast<TTree *>(fskim.Get("Events"));
    fskim.Close();

    // If skim has no Events tree, copy the whole file and return
    if (!skimEvents)
    {
        std::filesystem::copy_file(skimFile, outFriend, std::filesystem::copy_options::overwrite_existing);
        delete hClone;
        return true;
    }

    ROOT::RDataFrame df("Events", skimFile);
    auto df2 = df.Define("isGoodMuon_mva",
                         "(Muon_pt>15.0) && (abs(Muon_eta)<2.4) && (abs(Muon_dxy)<0.05) && (abs(Muon_dz)<0.1) && (abs(Muon_sip3d)<8) && Muon_mediumId && (Muon_miniPFRelIso_all<0.4) && (Muon_jetDF<0.2480) && (Muon_promptMVA > 0.64)")
                    .Define("nGoodMuon_mva", "ROOT::VecOps::Sum(isGoodMuon_mva)")
                    .Define("leadingMuonIdx", "Nonzero(isGoodMuon_mva)[0]")
                    .Define("subleadingMuonIdx", "Nonzero(isGoodMuon_mva)[1]");

    df2.Snapshot("Events", outFriend);

    // copy genWeightSum histogram into output file
    {
        TFile fout(outFriend.c_str(), "UPDATE");
        if (!fout.IsOpen())
        {
            rdfWS_utility::messageWARN("completeBranches_oneFile", "Failed to open output file for genWeightSum: " + outFriend);
            delete hClone;
            return false;
        }
        hClone->SetDirectory(&fout);
        hClone->Write("genWeightSum", TObject::kOverwrite);
    }
    delete hClone;
    return true;
}

int main(int argc, char *argv[])
{
    std::string skimmedFile, outFile;
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
                                       "Unknown argument: " + arg + ". Expected: --skimmed --outFile");
            return 1;
        }
    }

    if (skimmedFile.empty() || outFile.empty())
    {
        rdfWS_utility::messageERROR(
            "completeBranches_oneFile",
            "Input template wrong. Example: \"completeBranches_oneFile --skimmed <skimmed file> --outFile <output file>\"");
        return 1;
    }

    if (!process(skimmedFile, outFile))
        return 1;
    return 0;
}
