#include <functional>
#include <stdexcept>
#include <string>
#include <map>
#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

// for JES
std::shared_ptr<const correction::CorrectionSet> JEC_Data_cset;
std::shared_ptr<const correction::Correction> JEC_Data_corrL2;
std::shared_ptr<const correction::Correction> JEC_Data_corrL2L3;

std::shared_ptr<const correction::Correction> getDataCorrectionOrThrow(const std::vector<std::string>& names)
{
    std::string tried;
    for (const auto& name : names)
    {
        if (!tried.empty()) tried += ", ";
        tried += name;
        try {
            return JEC_Data_cset->at(name);
        } catch (const std::exception&) {
        }
    }
    throw std::runtime_error("JEC_Data correction not found. Tried: " + tried);
}

void JEC_Data_init(const std::string& era) { 
    std::map<std::string, std::string> JESFiles = {
        {"Run2024C", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/jet_jerc.json"},
        {"Run2024D", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/jet_jerc.json"},
        {"Run2024E", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/jet_jerc.json"},
        {"Run2024F", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/jet_jerc.json"},
        {"Run2024G", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/jet_jerc.json"},
        {"Run2024H", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/jet_jerc.json"},
        {"Run2024I", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/jet_jerc.json"},

        {"Run2022C", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-22CDSep23-Summer22-NanoAODv12/jet_jerc.json"},
        {"Run2022D", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-22CDSep23-Summer22-NanoAODv12/jet_jerc.json"},
        {"Run2022E", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/jet_jerc.json"},
        {"Run2022F", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/jet_jerc.json"},
        {"Run2022G", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/jet_jerc.json"},
        {"Run2023C", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-23CSep23-Summer23-NanoAODv12/jet_jerc.json"},
        {"Run2023D", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-23DSep23-Summer23BPix-NanoAODv12/jet_jerc.json"},
    };
    auto fileIt = JESFiles.find(era);
    if (fileIt == JESFiles.end()) {
        throw std::runtime_error("Unknown era for JEC data file: " + era);
    }
    JEC_Data_cset = correction::CorrectionSet::from_file(fileIt->second);

    static const std::map<std::string, std::string> L2Tab{
        {"Run2024C", "Summer24Prompt24_V2_DATA_L2Relative_AK4PFPuppi"},
        {"Run2024D", "Summer24Prompt24_V2_DATA_L2Relative_AK4PFPuppi"}, 
        {"Run2024E", "Summer24Prompt24_V2_DATA_L2Relative_AK4PFPuppi"},
        {"Run2024F", "Summer24Prompt24_V2_DATA_L2Relative_AK4PFPuppi"},
        {"Run2024G", "Summer24Prompt24_V2_DATA_L2Relative_AK4PFPuppi"},
        {"Run2024H", "Summer24Prompt24_V2_DATA_L2Relative_AK4PFPuppi"},
        {"Run2024I", "Summer24Prompt24_V2_DATA_L2Relative_AK4PFPuppi"},

        {"Run2022C", "Summer22_22Sep2023_V4_DATA_L2Relative_AK4PFPuppi"},
        {"Run2022D", "Summer22_22Sep2023_V4_DATA_L2Relative_AK4PFPuppi"},
        {"Run2022E", "Summer22EE_22Sep2023_V4_DATA_L2Relative_AK4PFPuppi"},
        {"Run2022F", "Summer22EE_22Sep2023_V4_DATA_L2Relative_AK4PFPuppi"},
        {"Run2022G", "Summer22EE_22Sep2023_V4_DATA_L2Relative_AK4PFPuppi"},
        {"Run2023C", "Summer23Prompt23_V4_DATA_L2Relative_AK4PFPuppi"},
        {"Run2023D", "Summer23BPixPrompt23_V4_DATA_L2Relative_AK4PFPuppi"},
    };
    auto it = L2Tab.find(era);
    if (it == L2Tab.end()) {
        throw std::runtime_error("Unknown era for L2: " + era);
    }
    JEC_Data_corrL2 = getDataCorrectionOrThrow({it->second});

    std::map<std::string, std::string> L2L3Tab{
        {"Run2024C", "Summer24Prompt24_V2_DATA_L2L3Residual_AK4PFPuppi"},
        {"Run2024D", "Summer24Prompt24_V2_DATA_L2L3Residual_AK4PFPuppi"},
        {"Run2024E", "Summer24Prompt24_V2_DATA_L2L3Residual_AK4PFPuppi"},
        {"Run2024F", "Summer24Prompt24_V2_DATA_L2L3Residual_AK4PFPuppi"},
        {"Run2024G", "Summer24Prompt24_V2_DATA_L2L3Residual_AK4PFPuppi"},
        {"Run2024H", "Summer24Prompt24_V2_DATA_L2L3Residual_AK4PFPuppi"},
        {"Run2024I", "Summer24Prompt24_V2_DATA_L2L3Residual_AK4PFPuppi"},

        {"Run2022C", "Summer22_22Sep2023_V4_DATA_L2L3Residual_AK4PFPuppi"},
        {"Run2022D", "Summer22_22Sep2023_V4_DATA_L2L3Residual_AK4PFPuppi"},
        {"Run2022E", "Summer22EE_22Sep2023_V4_DATA_L2L3Residual_AK4PFPuppi"},
        {"Run2022F", "Summer22EE_22Sep2023_V4_DATA_L2L3Residual_AK4PFPuppi"},
        {"Run2022G", "Summer22EE_22Sep2023_V4_DATA_L2L3Residual_AK4PFPuppi"},
        {"Run2023C", "Summer23Prompt23_V4_DATA_L2L3Residual_AK4PFPuppi"},
        {"Run2023D", "Summer23BPixPrompt23_V4_DATA_L2L3Residual_AK4PFPuppi"},
    };
    it = L2L3Tab.find(era);
    if (it == L2L3Tab.end()) {
        throw std::runtime_error("Unknown era for resolution: " + era);
    }
    JEC_Data_corrL2L3 = getDataCorrectionOrThrow({it->second});
}

// syst: "central" for nominal input
// flavour: b=5, c=4, udsg=0
ROOT::VecOps::RVec<float> get_JES_corr_pt(const ROOT::VecOps::RVec<float>& Var_to_correct, unsigned int run, const ROOT::VecOps::RVec<float>& Jet_pt, const ROOT::VecOps::RVec<float>& Jet_rawFactor, const ROOT::VecOps::RVec<float>& Jet_eta, const ROOT::VecOps::RVec<float>& Jet_phi)
{
    int nJet = Jet_pt.size();
    ROOT::VecOps::RVec<float> corrected(nJet);
    for (int i=0; i < nJet; i++)
    {
        float corrVar = Var_to_correct[i] * (1-Jet_rawFactor[i]);
        float corrPt = Jet_pt[i] * (1-Jet_rawFactor[i]);
        auto cL2 = JEC_Data_corrL2->evaluate({Jet_eta[i], Jet_phi[i], corrPt});
        corrVar *= cL2;
        corrPt *= cL2;
        auto cL2L3 = JEC_Data_corrL2L3->evaluate({float(run), Jet_eta[i], corrPt});
        corrVar *= cL2L3;
        corrected[i] = corrVar;
    }
    return corrected;
}

ROOT::VecOps::RVec<float> get_JES_corr_pt_bpix(const ROOT::VecOps::RVec<float>& Var_to_correct, unsigned int run, const ROOT::VecOps::RVec<float>& Jet_pt, const ROOT::VecOps::RVec<float>& Jet_rawFactor, const ROOT::VecOps::RVec<float>& Jet_eta, const ROOT::VecOps::RVec<float>& Jet_phi)
{
    int nJet = Jet_pt.size();
    ROOT::VecOps::RVec<float> corrected(nJet);
    for (int i=0; i < nJet; i++)
    {
        float corrVar = Var_to_correct[i] * (1-Jet_rawFactor[i]);
        float corrPt = Jet_pt[i] * (1-Jet_rawFactor[i]);
        auto cL2 = JEC_Data_corrL2->evaluate({Jet_eta[i], Jet_phi[i], corrPt});
        corrVar *= cL2;
        corrPt *= cL2;
        auto cL2L3 = JEC_Data_corrL2L3->evaluate({float(run), Jet_eta[i], corrPt});
        corrVar *= cL2L3;
        corrected[i] = corrVar;
    }
    return corrected;
}

ROOT::VecOps::RVec<float> get_JES_corr_pt_2023C(const ROOT::VecOps::RVec<float>& Var_to_correct, unsigned int run, const ROOT::VecOps::RVec<float>& Jet_pt, const ROOT::VecOps::RVec<float>& Jet_rawFactor, const ROOT::VecOps::RVec<float>& Jet_eta, const ROOT::VecOps::RVec<float>& Jet_phi)
{
    int nJet = Jet_pt.size();
    ROOT::VecOps::RVec<float> corrected(nJet);
    for (int i=0; i < nJet; i++)
    {
        float corrVar = Var_to_correct[i] * (1-Jet_rawFactor[i]);
        float corrPt = Jet_pt[i] * (1-Jet_rawFactor[i]);
        auto cL2 = JEC_Data_corrL2->evaluate({Jet_eta[i], corrPt});
        corrVar *= cL2;
        corrPt *= cL2;
        auto cL2L3 = JEC_Data_corrL2L3->evaluate({float(run), Jet_eta[i], corrPt});
        corrVar *= cL2L3;
        corrected[i] = corrVar;
    }
    return corrected;
}

ROOT::VecOps::RVec<float> get_JES_corr_pt_v12(const ROOT::VecOps::RVec<float>& Var_to_correct, unsigned int run, const ROOT::VecOps::RVec<float>& Jet_pt, const ROOT::VecOps::RVec<float>& Jet_rawFactor, const ROOT::VecOps::RVec<float>& Jet_eta, const ROOT::VecOps::RVec<float>& Jet_phi)
{
    int nJet = Jet_pt.size();
    ROOT::VecOps::RVec<float> corrected(nJet);
    for (int i=0; i < nJet; i++)
    {
        float corrVar = Var_to_correct[i] * (1-Jet_rawFactor[i]);
        float corrPt = Jet_pt[i] * (1-Jet_rawFactor[i]);
        auto cL2 = JEC_Data_corrL2->evaluate({Jet_eta[i], corrPt});
        corrVar *= cL2;
        corrPt *= cL2;
        auto cL2L3 = JEC_Data_corrL2L3->evaluate({float(run), Jet_eta[i], corrPt});
        corrVar *= cL2L3;
        corrected[i] = corrVar;
    }
    return corrected;
}
