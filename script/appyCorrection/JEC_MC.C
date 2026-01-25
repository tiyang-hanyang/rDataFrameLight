#include <functional>
#include <stdexcept>
#include <string>
#include <map>
#include <algorithm>
#include <cmath>
#include <memory>

// for JES
std::shared_ptr<const correction::CorrectionSet> JEC_MC_cset;
std::shared_ptr<const correction::Correction> JEC_MC_corrL2;

// for JER
std::shared_ptr<const correction::Correction> JEC_MC_corrReso;
std::shared_ptr<const correction::Correction> JEC_MC_corrScale;
std::shared_ptr<const correction::CorrectionSet> JEC_MC_smearSet;
std::shared_ptr<const correction::Correction> JEC_MC_smearMenu;

void JEC_MC_init(const std::string& era) { 
    std::map<std::string, std::string> JESFiles = {
        {"RunIII2024Summer24NanoAODv15", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/2024_Summer24/jet_jerc.json"},
        {"Run3Summer23BPixNanoAODv12", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/2023_Summer23BPix/jet_jerc.json"},

        {"Run3Summer22NanoAODv12", "/home/tiyang/public/rDataFrameLight_git/correction/JME/Run3-22CDSep23-Summer22-NanoAODv12/jet_jerc.json"},
        {"Run3Summer22EENanoAODv12", "/home/tiyang/public/rDataFrameLight_git/correction/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/jet_jerc.json"},
        {"Run3Summer23NanoAODv12", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/2023_Summer23/jet_jerc.json"},
    };
    JEC_MC_cset = correction::CorrectionSet::from_file(JESFiles[era]);

    // for JES
    static const std::map<std::string, std::string> L2Tab{
        {"RunIII2024Summer24NanoAODv15", "Summer24Prompt24_V1_MC_L2Relative_AK4PFPuppi"},
        {"Run3Summer23BPixNanoAODv12", "Summer23BPixPrompt23_V3_MC_L2Relative_AK4PFPuppi"},

        {"Run3Summer22NanoAODv12", "Summer22_22Sep2023_V3_MC_L2Relative_AK4PFPuppi"},
        {"Run3Summer22EENanoAODv12", "Summer22EE_22Sep2023_V3_MC_L2Relative_AK4PFPuppi"},
        {"Run3Summer23NanoAODv12", "Summer23Prompt23_V2_MC_L2Relative_AK4PFPuppi"},
    };
    auto it = L2Tab.find(era);
    if (it == L2Tab.end()) {
        throw std::runtime_error("Unknown era for L2: " + era);
    }
    JEC_MC_corrL2 = JEC_MC_cset->at(it->second);


    // for JER nom
    static const std::map<std::string, std::string> ResoTab{
        {"RunIII2024Summer24NanoAODv15", "Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi"},
        {"Run3Summer23BPixNanoAODv12", "Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi"},

        {"Run3Summer23NanoAODv12", "Summer23Prompt23_RunCv1234_JRV1_MC_PtResolution_AK4PFPuppi"},
        {"Run3Summer22EENanoAODv12", "Summer22EE_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi"},
        {"Run3Summer22NanoAODv12", "Summer22_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi"},
    };
    static const std::map<std::string, std::string> ScaleTab{
        {"RunIII2024Summer24NanoAODv15", "Summer23BPixPrompt23_RunD_JRV1_MC_ScaleFactor_AK4PFPuppi"},
        {"Run3Summer23BPixNanoAODv12", "Summer23BPixPrompt23_RunD_JRV1_MC_ScaleFactor_AK4PFPuppi"},

        {"Run3Summer23NanoAODv12", "Summer23Prompt23_RunCv1234_JRV1_MC_ScaleFactor_AK4PFPuppi"},
        {"Run3Summer22EENanoAODv12", "Summer22EE_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi"},
        {"Run3Summer22NanoAODv12", "Summer22_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi"},
    };
    it = ResoTab.find(era);
    if (it == ResoTab.end()) {
        throw std::runtime_error("Unknown era for resolution: " + era);
    }
    JEC_MC_corrReso = JEC_MC_cset->at(it->second);
    it = ScaleTab.find(era);
    if (it == ScaleTab.end()) {
        throw std::runtime_error("Unknown era for resolution: " + era);
    }
    JEC_MC_corrScale = JEC_MC_cset->at(it->second);

    // smear menu
    std::string smearSetFile = "/home/tiyang/public/rDataFrameLight_git/source/script/jer_smear.json";
    JEC_MC_smearSet = correction::CorrectionSet::from_file(smearSetFile);
    JEC_MC_smearMenu = JEC_MC_smearSet->at("JERSmear");
}

// syst: "central" for nominal input
// flavour: b=5, c=4, udsg=0
ROOT::VecOps::RVec<float> get_JES_corr_pt(const ROOT::VecOps::RVec<float>& Jet_pt, const ROOT::VecOps::RVec<float>& Jet_rawFactor, const ROOT::VecOps::RVec<float>& Jet_eta, const ROOT::VecOps::RVec<float>& Jet_phi)
{
    int nJet = Jet_pt.size();
    ROOT::VecOps::RVec<float> rawPt(nJet);
    for (int i=0; i < nJet; i++)
    {
        float rawPtVal = Jet_pt[i] * (1-Jet_rawFactor[i]);
        auto cL2 = JEC_MC_corrL2->evaluate({Jet_eta[i], Jet_phi[i], rawPtVal});
        rawPtVal *= cL2;
        rawPt[i] = rawPtVal;
    }
    return rawPt;
}

ROOT::VecOps::RVec<float> get_JES_corr_pt_v12(const ROOT::VecOps::RVec<float>& Jet_pt, const ROOT::VecOps::RVec<float>& Jet_rawFactor, const ROOT::VecOps::RVec<float>& Jet_eta, const ROOT::VecOps::RVec<float>& Jet_phi)
{
    int nJet = Jet_pt.size();
    ROOT::VecOps::RVec<float> rawPt(nJet);
    for (int i=0; i < nJet; i++)
    {
        float rawPtVal = Jet_pt[i] * (1-Jet_rawFactor[i]);
        auto cL2 = JEC_MC_corrL2->evaluate({Jet_eta[i], rawPtVal});
        rawPtVal *= cL2;
        rawPt[i] = rawPtVal;
    }
    return rawPt;
}

ROOT::VecOps::RVec<float> get_JER_corr_pt(const ROOT::VecOps::RVec<float>& Jet_pt_JES, const ROOT::VecOps::RVec<float>& Jet_eta, const ROOT::VecOps::RVec<float>& Jet_phi,const ROOT::VecOps::RVec<short>& Jet_genJetIdx, const ROOT::VecOps::RVec<float>& GenJet_pt, const ROOT::VecOps::RVec<float>& GenJet_eta, const ROOT::VecOps::RVec<float>& GenJet_phi, float Rho, UChar_t event)
{
    int nJet = Jet_pt_JES.size();
    ROOT::VecOps::RVec<float> smearPt(nJet);
    for (int i=0; i < nJet; i++)
    {
        float origPtVal = Jet_pt_JES[i];
        float eta = Jet_eta[i];
        float phi = Jet_phi[i];

        // from the JERC correctionlib
        float reso = JEC_MC_corrReso->evaluate({ eta, origPtVal, Rho });
        float sf = JEC_MC_corrScale->evaluate({ eta, origPtVal, "nom"});

        // check genjet matching
        float genPtForSmear = -1.0;

        int genIdx = Jet_genJetIdx[i];
        if (genIdx >= 0)
        {
            float genPt = GenJet_pt[genIdx];
            float genEta = GenJet_eta[genIdx];
            float genPhi = GenJet_phi[genIdx];
            
            float dPhi = std::abs(genPhi - phi); 
            if (dPhi > 3.1415926) dPhi = 2*3.1415926 - dPhi;
            float dEta = genEta - eta;
            float dR2 = dPhi*dPhi + dEta*dEta;
            float dPt = std::abs(genPt - origPtVal);
            // auto dR = std::sqrt(dPhi*dPhi + dEta * dEta);
            if (dR2 < 0.04 && (dPt < 3.0*reso*origPtVal)) 
                genPtForSmear = genPt;
        }

        auto JERcorr = JEC_MC_smearMenu->evaluate({origPtVal, eta, genPtForSmear, Rho, event, reso, sf});
        smearPt[i] = origPtVal * JERcorr;
    }
    return smearPt;
}