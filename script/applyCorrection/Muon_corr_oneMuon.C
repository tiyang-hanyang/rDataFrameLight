#include <functional>
#include <stdexcept>
#include <string>
#include <map>
#include <memory>
#include <algorithm>

std::shared_ptr<const correction::CorrectionSet> Muon_oneMuon_id_cset;
std::shared_ptr<const correction::CorrectionSet> Muon_oneMuon_mva_cset;
std::shared_ptr<const correction::Correction> Muon_oneMuon_id_corr;
std::shared_ptr<const correction::Correction> Muon_oneMuon_mva_corr;
bool Muon_oneMuon_use_combined_single_sf = false;

namespace {
bool UsesPrivateSingleMuonSF_oneMuon(const std::string& era) {
    return era == "Run3Summer22NanoAODv12" ||
           era == "Run3Summer22EENanoAODv12" ||
           era == "Run3Summer23NanoAODv12" ||
           era == "Run3Summer23BPixNanoAODv12";
}

float ClampMuonEta_oneMuon(float eta) {
    return std::clamp(eta, -2.39f, 2.39f);
}

float ClampMuonPt_oneMuon(float pt) {
    const float maxPt = Muon_oneMuon_use_combined_single_sf ? 499.9f : 1000.0f;
    return std::clamp(pt, 15.1f, maxPt);
}
}

void Muon_corr_oneMuon_init(const std::string& era, const std::string& moduleDir = "") {
    Muon_oneMuon_use_combined_single_sf = UsesPrivateSingleMuonSF_oneMuon(era);

    const auto requireFile = [&](const std::map<std::string, std::string> &fileMap, const std::string &label) -> std::string {
        auto it = fileMap.find(era);
        if (it == fileMap.end()) {
            throw std::runtime_error("Unknown era for " + label + ": " + era);
        }
        return it->second;
    };

    std::string idFile;
    std::string mvaFile;

    if (Muon_oneMuon_use_combined_single_sf) {
        static const std::map<std::string, std::string> PrivateMuonFiles = {
            {"Run3Summer22NanoAODv12", "Muon_sf_WZ/muon_mvaTTH_2022_SF.json"},
            {"Run3Summer22EENanoAODv12", "Muon_sf_WZ/muon_mvaTTH_2022EE_SF.json"},
            {"Run3Summer23NanoAODv12", "Muon_sf_WZ/muon_mvaTTH_2023.json"},
            {"Run3Summer23BPixNanoAODv12", "Muon_sf_WZ/muon_mvaTTH_2023BPix.json"},
        };
        const std::string relativeFile = requireFile(PrivateMuonFiles, "private muon correction");
        if (moduleDir.empty()) {
            throw std::runtime_error("moduleDir is required for private muon correction: " + era);
        }
        idFile = moduleDir + "/" + relativeFile;
        mvaFile = idFile;
    } else {
        static const std::map<std::string, std::string> OfficialMuonFiles = {
            {"RunIII2024Summer24NanoAODv15", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/MUO/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/muon_Z.json"},
        };
        idFile = requireFile(OfficialMuonFiles, "official muon correction");
        mvaFile = idFile;
    }

    Muon_oneMuon_id_cset = correction::CorrectionSet::from_file(idFile);
    Muon_oneMuon_mva_cset = correction::CorrectionSet::from_file(mvaFile);

    if (Muon_oneMuon_use_combined_single_sf) {
        static const std::string combinedName = "NUM_TightmvaTTH_DEN_LooseMuons";
        Muon_oneMuon_id_corr = Muon_oneMuon_id_cset->at(combinedName);
        Muon_oneMuon_mva_corr = Muon_oneMuon_mva_cset->at(combinedName);
    } else {
        static const std::map<std::string, std::string> IDTab{
            {"RunIII2024Summer24NanoAODv15", "NUM_MediumID_DEN_TrackerMuons"},
        };
        auto it = IDTab.find(era);
        if (it == IDTab.end()) {
            throw std::runtime_error("Unknown era for muon ID SF: " + era);
        }
        Muon_oneMuon_id_corr = Muon_oneMuon_id_cset->at(it->second);

        static const std::map<std::string, std::string> MVATab{
            {"RunIII2024Summer24NanoAODv15", "NUM_promptMVA_WP64ID_DEN_MediumID"},
        };
        it = MVATab.find(era);
        if (it == MVATab.end()) {
            throw std::runtime_error("Unknown era for muon MVA SF: " + era);
        }
        Muon_oneMuon_mva_corr = Muon_oneMuon_mva_cset->at(it->second);
    }
}

ROOT::VecOps::RVec<float> MuonIDScale_oneMuon(const ROOT::VecOps::RVec<float>& eta, const ROOT::VecOps::RVec<float>& pt, const std::string& syst = "nominal")
{
    ROOT::VecOps::RVec<float> MuonIDSF(eta.size());
    for (size_t i = 0; i < eta.size(); ++i)
    {
        float safePt = ClampMuonPt_oneMuon(pt[i]);
        float safeEta = ClampMuonEta_oneMuon(eta[i]);
        MuonIDSF[i] = Muon_oneMuon_id_corr->evaluate(std::vector<correction::Variable::Type>{safeEta, safePt, syst});
    }
    return MuonIDSF;
}

ROOT::VecOps::RVec<float> MuonMVAScale_oneMuon(const ROOT::VecOps::RVec<float>& eta, const ROOT::VecOps::RVec<float>& pt, const std::string& syst = "nominal")
{
    ROOT::VecOps::RVec<float> MuonMVASF(eta.size());
    for (size_t i = 0; i < eta.size(); ++i)
    {
        float safePt = ClampMuonPt_oneMuon(pt[i]);
        float safeEta = ClampMuonEta_oneMuon(eta[i]);
        MuonMVASF[i] = Muon_oneMuon_mva_corr->evaluate(std::vector<correction::Variable::Type>{safeEta, safePt, syst});
    }
    return MuonMVASF;
}

ROOT::VecOps::RVec<float> UnitMuonScale_oneMuon(const ROOT::VecOps::RVec<float>& eta)
{
    return ROOT::VecOps::RVec<float>(eta.size(), 1.0f);
}

float EventSingleMuonSF_oneMuon(const ROOT::VecOps::RVec<float>& SF,
                                unsigned int leadingMuonIdx)
{
    return SF[leadingMuonIdx];
}

float EventMuonSF_oneMuon(const ROOT::VecOps::RVec<float>& IDSF,
                          const ROOT::VecOps::RVec<float>& IsoSF,
                          const ROOT::VecOps::RVec<float>& MVASF,
                          unsigned int leadingMuonIdx)
{
    return IDSF[leadingMuonIdx] * IsoSF[leadingMuonIdx] * MVASF[leadingMuonIdx];
}
