#include <functional>
#include <stdexcept>
#include <string>
#include <map>
#include <algorithm>
#include <cmath>
#include <memory>

// for JES
std::shared_ptr<const correction::CorrectionSet> JEC_Data_cset;
std::shared_ptr<const correction::Correction> JEC_Data_corrL2;
std::shared_ptr<const correction::Correction> JEC_Data_corrL2L3;

void JEC_Data_init(const std::string& era) { 
    std::map<std::string, std::string> JESFiles = {
        {"Run2024C", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/2024_Summer24/jet_jerc.json"},
        {"Run2024D", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/2024_Summer24/jet_jerc.json"},
        {"Run2024E", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/2024_Summer24/jet_jerc.json"},
        {"Run2024F", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/2024_Summer24/jet_jerc.json"},
        {"Run2024G", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/2024_Summer24/jet_jerc.json"},
        {"Run2024H", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/2024_Summer24/jet_jerc.json"},
        {"Run2024I", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/2024_Summer24/jet_jerc.json"},

        {"Run2022C", "/home/tiyang/public/rDataFrameLight_git/correction/JME/Run3-22CDSep23-Summer22-NanoAODv12/jet_jerc.json"},
        {"Run2022D", "/home/tiyang/public/rDataFrameLight_git/correction/JME/Run3-22CDSep23-Summer22-NanoAODv12/jet_jerc.json"},
        {"Run2022E", "/home/tiyang/public/rDataFrameLight_git/correction/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/jet_jerc.json"},
        {"Run2022F", "/home/tiyang/public/rDataFrameLight_git/correction/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/jet_jerc.json"},
        {"Run2022G", "/home/tiyang/public/rDataFrameLight_git/correction/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/jet_jerc.json"},
        {"Run2023C", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/2023_Summer23/jet_jerc.json"},
        {"Run2023D", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/2023_Summer23BPix/jet_jerc.json"},
    };
    JEC_Data_cset = correction::CorrectionSet::from_file(JESFiles[era]);

    static const std::map<std::string, std::string> L2Tab{
        {"Run2024C", "Summer24Prompt24_V1_DATA_L2Relative_AK4PFPuppi"},
        {"Run2024D", "Summer24Prompt24_V1_DATA_L2Relative_AK4PFPuppi"}, 
        {"Run2024E", "Summer24Prompt24_V1_DATA_L2Relative_AK4PFPuppi"},
        {"Run2024F", "Summer24Prompt24_V1_DATA_L2Relative_AK4PFPuppi"},
        {"Run2024G", "Summer24Prompt24_V1_DATA_L2Relative_AK4PFPuppi"},
        {"Run2024H", "Summer24Prompt24_V1_DATA_L2Relative_AK4PFPuppi"},
        {"Run2024I", "Summer24Prompt24_V1_DATA_L2Relative_AK4PFPuppi"},

        {"Run2022C", "Summer22_22Sep2023_RunCD_V3_DATA_L2Relative_AK4PFPuppi"},
        {"Run2022D", "Summer22_22Sep2023_RunCD_V3_DATA_L2Relative_AK4PFPuppi"},
        {"Run2022E", "Summer22EE_22Sep2023_RunE_V3_DATA_L2Relative_AK4PFPuppi"},
        {"Run2022F", "Summer22EE_22Sep2023_RunF_V3_DATA_L2Relative_AK4PFPuppi"},
        {"Run2022G", "Summer22EE_22Sep2023_RunG_V3_DATA_L2Relative_AK4PFPuppi"},
        {"Run2023C", "Summer23Prompt23_V2_DATA_L2Relative_AK4PFPuppi"},
        {"Run2023D", "Summer23BPixPrompt23_V3_DATA_L2Relative_AK4PFPuppi"},
    };
    auto it = L2Tab.find(era);
    if (it == L2Tab.end()) {
        throw std::runtime_error("Unknown era for L2: " + era);
    }
    JEC_Data_corrL2 = JEC_Data_cset->at(it->second);

    std::map<std::string, std::string> L2L3Tab{
        {"Run2024C", "Summer24Prompt24_V1_DATA_L2L3Residual_AK4PFPuppi"},
        {"Run2024D", "Summer24Prompt24_V1_DATA_L2L3Residual_AK4PFPuppi"},
        {"Run2024E", "Summer24Prompt24_V1_DATA_L2L3Residual_AK4PFPuppi"},
        {"Run2024F", "Summer24Prompt24_V1_DATA_L2L3Residual_AK4PFPuppi"},
        {"Run2024G", "Summer24Prompt24_V1_DATA_L2L3Residual_AK4PFPuppi"},
        {"Run2024H", "Summer24Prompt24_V1_DATA_L2L3Residual_AK4PFPuppi"},
        {"Run2024I", "Summer24Prompt24_V1_DATA_L2L3Residual_AK4PFPuppi"},

        {"Run2022C", "Summer22_22Sep2023_RunCD_V3_DATA_L2L3Residual_AK4PFPuppi"},
        {"Run2022D", "Summer22_22Sep2023_RunCD_V3_DATA_L2L3Residual_AK4PFPuppi"},
        {"Run2022E", "Summer22EE_22Sep2023_RunE_V3_DATA_L2L3Residual_AK4PFPuppi"},
        {"Run2022F", "Summer22EE_22Sep2023_RunF_V3_DATA_L2L3Residual_AK4PFPuppi"},
        {"Run2022G", "Summer22EE_22Sep2023_RunG_V3_DATA_L2L3Residual_AK4PFPuppi"},
        {"Run2023C", "Summer23Prompt23_V2_DATA_L2L3Residual_AK4PFPuppi"},
        {"Run2023D", "Summer23BPixPrompt23_V3_DATA_L2L3Residual_AK4PFPuppi"},
    };
    it = L2L3Tab.find(era);
    if (it == L2L3Tab.end()) {
        throw std::runtime_error("Unknown era for resolution: " + era);
    }
    JEC_Data_corrL2L3 = JEC_Data_cset->at(it->second);
}

// syst: "central" for nominal input
// flavour: b=5, c=4, udsg=0
ROOT::VecOps::RVec<float> get_JES_corr_pt(unsigned int run, const ROOT::VecOps::RVec<float>& Jet_pt, const ROOT::VecOps::RVec<float>& Jet_rawFactor, const ROOT::VecOps::RVec<float>& Jet_eta, const ROOT::VecOps::RVec<float>& Jet_phi)
{
    int nJet = Jet_pt.size();
    ROOT::VecOps::RVec<float> rawPt(nJet);
    for (int i=0; i < nJet; i++)
    {
        float rawPtVal = Jet_pt[i] * (1-Jet_rawFactor[i]);
        auto cL2 = JEC_Data_corrL2->evaluate({Jet_eta[i], Jet_phi[i], rawPtVal});
        rawPtVal *= cL2;
        auto cL2L3 = JEC_Data_corrL2L3->evaluate({float(run), Jet_eta[i], rawPtVal});
        rawPtVal *= cL2L3;
        rawPt[i] = rawPtVal;
    }
    return rawPt;
}

ROOT::VecOps::RVec<float> get_JES_corr_pt_2023C(unsigned int run, const ROOT::VecOps::RVec<float>& Jet_pt, const ROOT::VecOps::RVec<float>& Jet_rawFactor, const ROOT::VecOps::RVec<float>& Jet_eta, const ROOT::VecOps::RVec<float>& Jet_phi)
{
    int nJet = Jet_pt.size();
    ROOT::VecOps::RVec<float> rawPt(nJet);
    for (int i=0; i < nJet; i++)
    {
        float rawPtVal = Jet_pt[i] * (1-Jet_rawFactor[i]);
        auto cL2 = JEC_Data_corrL2->evaluate({Jet_eta[i], rawPtVal});
        rawPtVal *= cL2;
        auto cL2L3 = JEC_Data_corrL2L3->evaluate({float(run), Jet_eta[i], rawPtVal});
        rawPtVal *= cL2L3;
        rawPt[i] = rawPtVal;
    }
    return rawPt;
}

ROOT::VecOps::RVec<float> get_JES_corr_pt_v12(unsigned int run, const ROOT::VecOps::RVec<float>& Jet_pt, const ROOT::VecOps::RVec<float>& Jet_rawFactor, const ROOT::VecOps::RVec<float>& Jet_eta, const ROOT::VecOps::RVec<float>& Jet_phi)
{
    int nJet = Jet_pt.size();
    ROOT::VecOps::RVec<float> rawPt(nJet);
    for (int i=0; i < nJet; i++)
    {
        float rawPtVal = Jet_pt[i] * (1-Jet_rawFactor[i]);
        auto cL2 = JEC_Data_corrL2->evaluate({Jet_eta[i], rawPtVal});
        rawPtVal *= cL2;
        auto cL2L3 = JEC_Data_corrL2L3->evaluate({Jet_eta[i], rawPtVal});
        rawPtVal *= cL2L3;
        rawPt[i] = rawPtVal;
    }
    return rawPt;
}