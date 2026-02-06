#include <functional>
#include <stdexcept>
#include <string>
#include <map>
#include <memory>
#include <algorithm>

std::shared_ptr<const correction::CorrectionSet> Muon_cset;
std::shared_ptr<const correction::Correction> Muon_id_corr;
// note that there is no any muon looseMiniIso_DEN_TightID, but only loose MiniIso_DEN_MediumID
std::shared_ptr<const correction::Correction> Muon_iso_corr;
std::shared_ptr<const correction::Correction> Muon_mva_corr;

void Muon_corr_init(const std::string& era) { 
    std::map<std::string, std::string> MuonFiles = {
        {"RunIII2024Summer24NanoAODv15", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/MUO/2024_Summer24/muon_Z.json"},
        {"Run3Summer23NanoAODv12", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/MUO/2023_Summer23BPix/muon_Z.json"},
    };
    Muon_cset = correction::CorrectionSet::from_file(MuonFiles[era]);

    // ID
    static const std::map<std::string, std::string> IDTab{
        {"RunIII2024Summer24NanoAODv15", "NUM_MediumID_DEN_TrackerMuons"},
        {"Run3Summer23NanoAODv12", "NUM_MediumID_DEN_TrackerMuons"},
    };
    auto it = IDTab.find(era);
    if (it == IDTab.end()) {
        throw std::runtime_error("Unknown era for L2: " + era);
    }
    Muon_id_corr = Muon_cset->at(it->second);

    // ISO
    static const std::map<std::string, std::string> ISOTab{
        {"RunIII2024Summer24NanoAODv15", "NUM_LooseMiniIso_DEN_MediumID"},
        {"Run3Summer23NanoAODv12", "NUM_LooseMiniIso_DEN_MediumID"},
    };
    it = ISOTab.find(era);
    if (it == ISOTab.end()) {
        throw std::runtime_error("Unknown era for L2: " + era);
    }
    Muon_iso_corr = Muon_cset->at(it->second);

    // MVA
    static const std::map<std::string, std::string> MVATab{
        {"RunIII2024Summer24NanoAODv15", "NUM_promptMVA_WP64ID_DEN_MediumID"},
        {"Run3Summer23NanoAODv12", "NUM_promptMVA_WP64ID_DEN_MediumID"},
    };
    it = MVATab.find(era);
    if (it == MVATab.end()) {
        throw std::runtime_error("Unknown era for L2: " + era);
    }
    Muon_mva_corr = Muon_cset->at(it->second);

}

ROOT::VecOps::RVec<float> MuonIDScale(const ROOT::VecOps::RVec<float>& eta, const ROOT::VecOps::RVec<float>& pt)
{
    ROOT::VecOps::RVec<float> MuonIDSF(eta.size());
    for (size_t i = 0; i < eta.size(); ++i)
    {
        float safePt = std::clamp(pt[i], 15.1f, 1000.0f);
        float safeEta = std::clamp(eta[i], -2.39f, 2.39f);
        // float safe_eta = std::clamp(eta[i], -2.4f, 2.4f);
        MuonIDSF[i] = Muon_id_corr->evaluate(std::vector<correction::Variable::Type>{safeEta, safePt, "nominal"});
    }
    return MuonIDSF;
}

ROOT::VecOps::RVec<float> MuonIsoScale(const ROOT::VecOps::RVec<float>& eta, const ROOT::VecOps::RVec<float>& pt)
{
    ROOT::VecOps::RVec<float> MuonIsoSF(eta.size());
    for (size_t i = 0; i < eta.size(); ++i)
    {
        float safePt = std::clamp(pt[i], 15.1f, 1000.0f);
        float safeEta = std::clamp(eta[i], -2.39f, 2.39f);
        // float safe_eta = std::clamp(eta[i], -2.4f, 2.4f);
        MuonIsoSF[i] = Muon_iso_corr->evaluate(std::vector<correction::Variable::Type>{safeEta, safePt, "nominal"});
    }
    return MuonIsoSF;
}

ROOT::VecOps::RVec<float> MuonMVAScale(const ROOT::VecOps::RVec<float>& eta, const ROOT::VecOps::RVec<float>& pt)
{
    ROOT::VecOps::RVec<float> MuonMVASF(eta.size());
    for (size_t i = 0; i < eta.size(); ++i)
    {
        float safePt = std::clamp(pt[i], 15.1f, 1000.0f);
        float safeEta = std::clamp(eta[i], -2.39f, 2.39f);
        // float safe_eta = std::clamp(eta[i], -2.4f, 2.4f);
        MuonMVASF[i] = Muon_mva_corr->evaluate(std::vector<correction::Variable::Type>{safeEta, safePt, "nominal"});
    }
    return MuonMVASF;
}


float EventMuonSF(const ROOT::VecOps::RVec<float>& IDSF,
                  const ROOT::VecOps::RVec<float>& IsoSF,
                  const ROOT::VecOps::RVec<float>& MVASF,
                  unsigned int leadingMuonIdx, unsigned int subleadingMuonIdx)
{
    if (leadingMuonIdx < 0 || subleadingMuonIdx < 0) return 1.0;
    return IDSF[leadingMuonIdx] * IDSF[subleadingMuonIdx] * IsoSF[leadingMuonIdx] * IsoSF[subleadingMuonIdx] * MVASF[leadingMuonIdx] * MVASF[subleadingMuonIdx];
}