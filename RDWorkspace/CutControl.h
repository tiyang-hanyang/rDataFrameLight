#ifndef CUT_CONTROL_H
#define CUT_CONTROL_H

#include "external/json.hpp"

#include <vector>
#include <tuple>
#include <string>
#include <functional>
#include <optional>

#include "ROOT/RDataFrame.hxx"

struct DsInfo
{
    int isMC=0;
    double lumi=9.451; // default value Run2023D
    double XS=6221.3; // default value DY
};

typedef std::tuple<std::string, std::string, std::string> cut_string;

class CutControl
{
    private:
        // achieve from json config files
        // 3-tuples of (operation, name, expression)
        // example
        // ["define", "leadingPt", "GoodMuon_pt[0]"], define a ordinary float
        // ["cut", "dimuon", "(nGoodMuon==2)"], make the cut of requiring exactly 2 muons, the name is only for debugging
        // ["TLVdefine", "leadingP4", "leadingPt,leadingEta,leadingPhi,leadingM"], define TLorentzVector for further usage, splitting the steps for reducing pre-defined lambda Define, e.g. find diMuon properties
        std::vector<cut_string> _steps; 

        // lambda function for applying the cut, constructed from json and store here.
        // lazy execution, only when the first time invoking applyCut, the _applyLambda is initialized
        mutable std::optional<std::function<ROOT::RDF::RNode(ROOT::RDF::RNode)>> _applyLambda;

        // for handling TLorentzVector needed case
        ROOT::RDF::RNode _defineTLV(ROOT::RDF::RNode origName, std::string name); 

    public:
        CutControl()=default;
        CutControl(const std::string& jsonFileName);
        CutControl(nlohmann::json json);
        // copy constructor, only copy the _steps, but not _applyLambda
        CutControl(const CutControl& origControl);
        ~CutControl()=default;

        void loadProcedure(const std::string& jsonFileName);
        void loadProcedure(nlohmann::json json);

        // can manually add new steps into THIS object and clear the _applyLambda, handle with care !!!
        void addStep(const cut_string moreStep); 
        // return a new CutControl object with additional CutControl stick after this 
        CutControl extend(const CutControl& addCon) const; 

        // core function, if _applyLambda not initialized, it will be initialized according to _steps
        // will warn if _step is empty
        ROOT::RDF::RNode applyCut(ROOT::RDF::RNode origRDF);

        // for parsing captured variables for TLorentzVector components
        std::vector<std::string> extractTLVComp(const std::string& TLVComp) const;

        // for debugging
        void printSteps();
};

// give convenient addition operation to concatenate cuts!
CutControl operator+(const CutControl& con1, const CutControl& con2);

#endif