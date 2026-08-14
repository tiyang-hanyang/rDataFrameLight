#include <functional>
#include <stdexcept>
#include <string>
#include <map>
#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

// for JES
std::shared_ptr<const correction::CorrectionSet> JEC_MC_cset;
std::shared_ptr<const correction::Correction> JEC_MC_corrL2;
std::map<std::string, std::shared_ptr<const correction::Correction>> JEC_MC_jesUncCorr;

// for JER
std::shared_ptr<const correction::Correction> JEC_MC_corrReso;
std::shared_ptr<const correction::Correction> JEC_MC_corrScale;
std::shared_ptr<const correction::CorrectionSet> JEC_MC_smearSet;
std::shared_ptr<const correction::Correction> JEC_MC_smearMenu;

int eventIdForJERSmear(ULong64_t event)
{
    return static_cast<int>(event & 0x7fffffffULL);
}

std::shared_ptr<const correction::Correction> getCorrectionOrThrow(const std::string& name)
{
    try {
        return JEC_MC_cset->at(name);
    } catch (const std::exception& e) {
        throw std::runtime_error("JEC_MC correction not found: " + name + " (" + e.what() + ")");
    }
}

void loadJESUncertainties(const std::string& prefix, const std::string& year)
{
    static const std::vector<std::string> commonSources{
        "FlavorQCD",
        "RelativeBal",
        "HF",
        "BBEC1",
        "EC2",
        "Absolute",
    };
    static const std::vector<std::string> yearSources{
        "Absolute",
        "HF",
        "EC2",
        "RelativeSample",
        "BBEC1",
    };
    for (const auto& source : commonSources) {
        const std::string branchName = "CMS_scale_j_" + source;
        const std::string corrName = prefix + "_MC_Regrouped_" + source + "_AK4PFPuppi";
        JEC_MC_jesUncCorr[branchName] = getCorrectionOrThrow(corrName);
    }
    for (const auto& source : yearSources) {
        const std::string branchName = "CMS_scale_j_" + source + "_" + year;
        const std::string corrName = prefix + "_MC_Regrouped_" + source + "_" + year + "_AK4PFPuppi";
        JEC_MC_jesUncCorr[branchName] = getCorrectionOrThrow(corrName);
    }
}

void JEC_MC_init(const std::string& era) { 
    std::map<std::string, std::string> JESFiles = {
        {"RunIII2024Summer24NanoAODv15", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/jet_jerc.json"},
        {"RunIII2024Summer24NanoAODv15_SSCR", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/jet_jerc.json"},
        {"RunIII2024Summer24NanoAODv15_AR", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/jet_jerc.json"},
        {"Run3Summer23BPixNanoAODv12", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-23DSep23-Summer23BPix-NanoAODv12/jet_jerc.json"},

        {"Run3Summer22NanoAODv12", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-22CDSep23-Summer22-NanoAODv12/jet_jerc.json"},
        {"Run3Summer22EENanoAODv12", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/jet_jerc.json"},
        {"Run3Summer23NanoAODv12", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-23CSep23-Summer23-NanoAODv12/jet_jerc.json"},
    };
    auto fileIt = JESFiles.find(era);
    if (fileIt == JESFiles.end()) {
        throw std::runtime_error("Unknown era for JEC/JER file: " + era);
    }
    JEC_MC_cset = correction::CorrectionSet::from_file(fileIt->second);

    // for JES
    static const std::map<std::string, std::string> L2Tab{
        {"RunIII2024Summer24NanoAODv15", "Summer24Prompt24_V2_MC_L2Relative_AK4PFPuppi"},
        {"RunIII2024Summer24NanoAODv15_SSCR", "Summer24Prompt24_V2_MC_L2Relative_AK4PFPuppi"},
        {"RunIII2024Summer24NanoAODv15_AR", "Summer24Prompt24_V2_MC_L2Relative_AK4PFPuppi"},
        {"Run3Summer23BPixNanoAODv12", "Summer23BPixPrompt23_V3_MC_L2Relative_AK4PFPuppi"},

        {"Run3Summer22NanoAODv12", "Summer22_22Sep2023_V3_MC_L2Relative_AK4PFPuppi"},
        {"Run3Summer22EENanoAODv12", "Summer22EE_22Sep2023_V3_MC_L2Relative_AK4PFPuppi"},
        {"Run3Summer23NanoAODv12", "Summer23Prompt23_V3_MC_L2Relative_AK4PFPuppi"},
    };
    auto it = L2Tab.find(era);
    if (it == L2Tab.end()) {
        throw std::runtime_error("Unknown era for L2: " + era);
    }
    JEC_MC_corrL2 = getCorrectionOrThrow(it->second);

    JEC_MC_jesUncCorr.clear();
    if (era == "RunIII2024Summer24NanoAODv15" || era == "RunIII2024Summer24NanoAODv15_SSCR" || era == "RunIII2024Summer24NanoAODv15_AR") {
        loadJESUncertainties("Summer24Prompt24_V2", "2024");
    } else if (era == "Run3Summer23BPixNanoAODv12") {
        loadJESUncertainties("Summer23BPixPrompt23_V3", "2023BPix");
    } else if (era == "Run3Summer23NanoAODv12") {
        loadJESUncertainties("Summer23Prompt23_V3", "2023");
    } else if (era == "Run3Summer22EENanoAODv12") {
        loadJESUncertainties("Summer22EE_22Sep2023_V3", "2022EE");
    } else if (era == "Run3Summer22NanoAODv12") {
        loadJESUncertainties("Summer22_22Sep2023_V3", "2022");
    }


    // for JER nom
    static const std::map<std::string, std::string> ResoTab{
        {"RunIII2024Summer24NanoAODv15", "Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi"},
        {"RunIII2024Summer24NanoAODv15_SSCR", "Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi"},
        {"RunIII2024Summer24NanoAODv15_AR", "Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi"},
        {"Run3Summer23BPixNanoAODv12", "Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi"},

        {"Run3Summer23NanoAODv12", "Summer23Prompt23_RunCv1234_JRV1_MC_PtResolution_AK4PFPuppi"},
        {"Run3Summer22EENanoAODv12", "Summer22EE_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi"},
        {"Run3Summer22NanoAODv12", "Summer22_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi"},
    };
    static const std::map<std::string, std::string> ScaleTab{
        {"RunIII2024Summer24NanoAODv15", "Summer23BPixPrompt23_RunD_JRV1_MC_ScaleFactor_AK4PFPuppi"},
        {"RunIII2024Summer24NanoAODv15_SSCR", "Summer23BPixPrompt23_RunD_JRV1_MC_ScaleFactor_AK4PFPuppi"},
        {"RunIII2024Summer24NanoAODv15_AR", "Summer23BPixPrompt23_RunD_JRV1_MC_ScaleFactor_AK4PFPuppi"},
        {"Run3Summer23BPixNanoAODv12", "Summer23BPixPrompt23_RunD_JRV1_MC_ScaleFactor_AK4PFPuppi"},

        {"Run3Summer23NanoAODv12", "Summer23Prompt23_RunCv1234_JRV1_MC_ScaleFactor_AK4PFPuppi"},
        {"Run3Summer22EENanoAODv12", "Summer22EE_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi"},
        {"Run3Summer22NanoAODv12", "Summer22_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi"},
    };
    it = ResoTab.find(era);
    if (it == ResoTab.end()) {
        throw std::runtime_error("Unknown era for resolution: " + era);
    }
    JEC_MC_corrReso = getCorrectionOrThrow(it->second);
    it = ScaleTab.find(era);
    if (it == ScaleTab.end()) {
        throw std::runtime_error("Unknown era for resolution: " + era);
    }
    JEC_MC_corrScale = getCorrectionOrThrow(it->second);

    // smear menu
    std::string smearSetFile = "/home/tiyang/public/rDataFrameLight_git/source/script/jer_smear.json";
    JEC_MC_smearSet = correction::CorrectionSet::from_file(smearSetFile);
    JEC_MC_smearMenu = JEC_MC_smearSet->at("JERSmear");
}

// syst: "central" for nominal input
// flavour: b=5, c=4, udsg=0
ROOT::VecOps::RVec<float> get_JES_corr_pt(const ROOT::VecOps::RVec<float>& Var_to_correct, const ROOT::VecOps::RVec<float>& Jet_pt, const ROOT::VecOps::RVec<float>& Jet_rawFactor, const ROOT::VecOps::RVec<float>& Jet_eta, const ROOT::VecOps::RVec<float>& Jet_phi)
{
    int nJet = Jet_pt.size();
    ROOT::VecOps::RVec<float> corrected(nJet);
    for (int i=0; i < nJet; i++)
    {
        float corrVar = Var_to_correct[i] * (1-Jet_rawFactor[i]);
        float corrPt = Jet_pt[i] * (1-Jet_rawFactor[i]);
        auto cL2 = JEC_MC_corrL2->evaluate({Jet_eta[i], Jet_phi[i], corrPt});
        corrVar *= cL2;
        corrected[i] = corrVar;
    }
    return corrected;
}

ROOT::VecOps::RVec<float> get_JES_corr_pt_v12(const ROOT::VecOps::RVec<float>& Var_to_correct, const ROOT::VecOps::RVec<float>& Jet_pt, const ROOT::VecOps::RVec<float>& Jet_rawFactor, const ROOT::VecOps::RVec<float>& Jet_eta, const ROOT::VecOps::RVec<float>& Jet_phi)
{
    int nJet = Jet_pt.size();
    ROOT::VecOps::RVec<float> corrected(nJet);
    for (int i=0; i < nJet; i++)
    {
        float corrVar = Var_to_correct[i] * (1-Jet_rawFactor[i]);
        float corrPt = Jet_pt[i] * (1-Jet_rawFactor[i]);
        auto cL2 = JEC_MC_corrL2->evaluate({Jet_eta[i], corrPt});
        corrVar *= cL2;
        corrected[i] = corrVar;
    }
    return corrected;
}

ROOT::VecOps::RVec<float> get_JES_uncertainty(const std::string& sourceName, const ROOT::VecOps::RVec<float>& Jet_eta, const ROOT::VecOps::RVec<float>& Jet_pt)
{
    auto corrIt = JEC_MC_jesUncCorr.find(sourceName);
    if (corrIt == JEC_MC_jesUncCorr.end()) {
        throw std::runtime_error("Unknown or unavailable JES uncertainty source: " + sourceName);
    }
    int nJet = Jet_pt.size();
    ROOT::VecOps::RVec<float> uncertainty(nJet);
    for (int i=0; i < nJet; i++)
    {
        uncertainty[i] = corrIt->second->evaluate({Jet_eta[i], Jet_pt[i]});
    }
    return uncertainty;
}

ROOT::VecOps::RVec<float> get_JER_corr(const ROOT::VecOps::RVec<float>& Jet_pt_JES, const ROOT::VecOps::RVec<float>& Jet_eta, const ROOT::VecOps::RVec<float>& Jet_phi,const ROOT::VecOps::RVec<short>& Jet_genJetIdx, const ROOT::VecOps::RVec<float>& GenJet_pt, const ROOT::VecOps::RVec<float>& GenJet_eta, const ROOT::VecOps::RVec<float>& GenJet_phi, float Rho, ULong64_t event, const std::string& syst = "nom")
{
    int nJet = Jet_pt_JES.size();
    ROOT::VecOps::RVec<float> smearCorr(nJet);
    for (int i=0; i < nJet; i++)
    {
        float origPtVal = Jet_pt_JES[i];
        float eta = Jet_eta[i];
        float phi = Jet_phi[i];

        // from the JERC correctionlib
        float reso = JEC_MC_corrReso->evaluate({ eta, origPtVal, Rho });
        float sf = JEC_MC_corrScale->evaluate({ eta, origPtVal, syst});

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

        auto JERcorr = JEC_MC_smearMenu->evaluate({origPtVal, eta, genPtForSmear, Rho, eventIdForJERSmear(event), reso, sf});
        smearCorr[i] = JERcorr;
    }
    return smearCorr;
}

ROOT::VecOps::RVec<float> get_JER_corr_pt(const ROOT::VecOps::RVec<float>& varToCorrect, const ROOT::VecOps::RVec<float>& Jet_pt_JES, const ROOT::VecOps::RVec<float>& Jet_eta, const ROOT::VecOps::RVec<float>& Jet_phi,const ROOT::VecOps::RVec<short>& Jet_genJetIdx, const ROOT::VecOps::RVec<float>& GenJet_pt, const ROOT::VecOps::RVec<float>& GenJet_eta, const ROOT::VecOps::RVec<float>& GenJet_phi, float Rho, ULong64_t event, const std::string& syst = "nom")
{
    return varToCorrect * get_JER_corr(Jet_pt_JES, Jet_eta, Jet_phi, Jet_genJetIdx, GenJet_pt, GenJet_eta, GenJet_phi, Rho, event, syst);
}
