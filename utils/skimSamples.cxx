#include "Utility.h"
#include "SkimControl.h"

#include "ROOT/RDataFrame.hxx"

int main(int argc, char *argv[])
{
    ROOT::EnableImplicitMT();
    rdfWS_utility::messageINFO("skimSamples", "Begin skimming jobs");

    // default skim config
    std::string skimConfig = "json/DiMuonSkim.json";
    if (argc > 1)
        skimConfig = argv[1];

    SkimControl mySkimmer(skimConfig);
    // can control what samples to skim
    // mySkimmer.turnOff(".*");
    // mySkimmer.turnOn("Muon0");
    mySkimmer.run();

    return 0;
}