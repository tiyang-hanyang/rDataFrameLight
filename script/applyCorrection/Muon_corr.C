#include <functional>
#include <stdexcept>
#include <string>
#include <map>
#include <memory>
#include <algorithm>

std::shared_ptr<const correction::CorrectionSet> Muon_id_cset;
std::shared_ptr<const correction::CorrectionSet> Muon_mva_cset;
std::shared_ptr<const correction::Correction> Muon_id_corr;
std::shared_ptr<const correction::Correction> Muon_mva_corr;
bool Muon_use_combined_single_sf = false;

namespace {
bool UsesPrivateSingleMuonSF(const std::string& era) {
    return era == "Run3Summer22NanoAODv12" ||
           era == "Run3Summer22EENanoAODv12" ||
           era == "Run3Summer23NanoAODv12" ||
           era == "Run3Summer23BPixNanoAODv12";
}

float ClampMuonEta(float eta) {
    return std::clamp(eta, -2.39f, 2.39f);
}

float ClampMuonPt(float pt) {
    // The private Run3 muon_mvaTTH jsons have a hard upper edge at 500 GeV with
    // flow="error", so stay just inside the last bin. The official 2024
    // correctionlib supports the broader range used previously.
    const float maxPt = Muon_use_combined_single_sf ? 499.9f : 1000.0f;
    return std::clamp(pt, 15.1f, maxPt);
}
}

void Muon_corr_init(const std::string& era, const std::string& moduleDir = "") {
    Muon_use_combined_single_sf = UsesPrivateSingleMuonSF(era);

    const auto requireFile = [&](const std::map<std::string, std::string> &fileMap, const std::string &label) -> std::string {
        auto it = fileMap.find(era);
        if (it == fileMap.end()) {
            throw std::runtime_error("Unknown era for " + label + ": " + era);
        }
        return it->second;
    };

    std::string idFile;
    std::string mvaFile;

    if (Muon_use_combined_single_sf) {
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
        // 2024 keeps using the official correctionlib; ISO is intentionally excluded
        // later when composing the event-level MuonScale.
        idFile = requireFile(OfficialMuonFiles, "official muon correction");
        mvaFile = idFile;
    }

    Muon_id_cset = correction::CorrectionSet::from_file(idFile);
    Muon_mva_cset = correction::CorrectionSet::from_file(mvaFile);

    if (Muon_use_combined_single_sf) {
        static const std::string combinedName = "NUM_TightmvaTTH_DEN_LooseMuons";
        Muon_id_corr = Muon_id_cset->at(combinedName);
        Muon_mva_corr = Muon_mva_cset->at(combinedName);
    } else {
        // ID
        static const std::map<std::string, std::string> IDTab{
            {"RunIII2024Summer24NanoAODv15", "NUM_MediumID_DEN_TrackerMuons"},
        };
        auto it = IDTab.find(era);
        if (it == IDTab.end()) {
            throw std::runtime_error("Unknown era for muon ID SF: " + era);
        }
        Muon_id_corr = Muon_id_cset->at(it->second);

        // MVA
        static const std::map<std::string, std::string> MVATab{
            {"RunIII2024Summer24NanoAODv15", "NUM_promptMVA_WP64ID_DEN_MediumID"},
        };
        it = MVATab.find(era);
        if (it == MVATab.end()) {
            throw std::runtime_error("Unknown era for muon MVA SF: " + era);
        }
        Muon_mva_corr = Muon_mva_cset->at(it->second);
    }

}

ROOT::VecOps::RVec<float> MuonIDScale(const ROOT::VecOps::RVec<float>& eta, const ROOT::VecOps::RVec<float>& pt, const std::string& syst = "nominal")
{
    ROOT::VecOps::RVec<float> MuonIDSF(eta.size());
    for (size_t i = 0; i < eta.size(); ++i)
    {
        float safePt = ClampMuonPt(pt[i]);
        float safeEta = ClampMuonEta(eta[i]);
        MuonIDSF[i] = Muon_id_corr->evaluate(std::vector<correction::Variable::Type>{safeEta, safePt, syst});
    }
    return MuonIDSF;
}

ROOT::VecOps::RVec<float> MuonIsoScale(const ROOT::VecOps::RVec<float>& eta, const ROOT::VecOps::RVec<float>& pt, const std::string& syst = "nominal")
{
    return ROOT::VecOps::RVec<float>(eta.size(), 1.0f);
}

ROOT::VecOps::RVec<float> MuonMVAScale(const ROOT::VecOps::RVec<float>& eta, const ROOT::VecOps::RVec<float>& pt, const std::string& syst = "nominal")
{
    ROOT::VecOps::RVec<float> MuonMVASF(eta.size());
    for (size_t i = 0; i < eta.size(); ++i)
    {
        float safePt = ClampMuonPt(pt[i]);
        float safeEta = ClampMuonEta(eta[i]);
        MuonMVASF[i] = Muon_mva_corr->evaluate(std::vector<correction::Variable::Type>{safeEta, safePt, syst});
    }
    return MuonMVASF;
}

ROOT::VecOps::RVec<float> UnitMuonScale(const ROOT::VecOps::RVec<float>& eta)
{
    return ROOT::VecOps::RVec<float>(eta.size(), 1.0f);
}

float EventSingleMuonSF(const ROOT::VecOps::RVec<float>& SF,
                        unsigned int leadingMuonIdx, unsigned int subleadingMuonIdx)
{
    if (leadingMuonIdx < 0 || subleadingMuonIdx < 0) return 1.0;
    return SF[leadingMuonIdx] * SF[subleadingMuonIdx];
}

float EventSingleMuonSFSelected(const ROOT::VecOps::RVec<float>& SF,
                                const ROOT::VecOps::RVec<int>& goodMuonMask,
                                const ROOT::VecOps::RVec<float>& pt,
                                float leadingPtThreshold = 20.0f,
                                float subleadingPtThreshold = 15.0f)
{
    float scale = 1.0f;
    float maxSelectedPt = -1.0f;
    bool hasSelectedMuon = false;
    for (size_t i = 0; i < SF.size() && i < goodMuonMask.size() && i < pt.size(); ++i)
    {
        if (!goodMuonMask[i] || pt[i] <= subleadingPtThreshold) continue;
        hasSelectedMuon = true;
        maxSelectedPt = std::max(maxSelectedPt, pt[i]);
        scale *= SF[i];
    }
    if (!hasSelectedMuon || maxSelectedPt <= leadingPtThreshold) return 1.0f;
    return scale;
}

float EventMuonSF(const ROOT::VecOps::RVec<float>& IDSF,
                  const ROOT::VecOps::RVec<float>& IsoSF,
                  const ROOT::VecOps::RVec<float>& MVASF,
                  unsigned int leadingMuonIdx, unsigned int subleadingMuonIdx)
{
    if (leadingMuonIdx < 0 || subleadingMuonIdx < 0) return 1.0;
    return IDSF[leadingMuonIdx] * IDSF[subleadingMuonIdx] * IsoSF[leadingMuonIdx] * IsoSF[subleadingMuonIdx] * MVASF[leadingMuonIdx] * MVASF[subleadingMuonIdx];
}

float EventMuonSFSelected(const ROOT::VecOps::RVec<float>& IDSF,
                          const ROOT::VecOps::RVec<float>& IsoSF,
                          const ROOT::VecOps::RVec<float>& MVASF,
                          const ROOT::VecOps::RVec<int>& goodMuonMask,
                          const ROOT::VecOps::RVec<float>& pt,
                          float leadingPtThreshold = 20.0f,
                          float subleadingPtThreshold = 15.0f)
{
    float scale = 1.0f;
    float maxSelectedPt = -1.0f;
    bool hasSelectedMuon = false;
    for (size_t i = 0;
         i < IDSF.size() && i < IsoSF.size() && i < MVASF.size() && i < goodMuonMask.size() && i < pt.size();
         ++i)
    {
        if (!goodMuonMask[i] || pt[i] <= subleadingPtThreshold) continue;
        hasSelectedMuon = true;
        maxSelectedPt = std::max(maxSelectedPt, pt[i]);
        scale *= IDSF[i] * IsoSF[i] * MVASF[i];
    }
    if (!hasSelectedMuon || maxSelectedPt <= leadingPtThreshold) return 1.0f;
    return scale;
}
