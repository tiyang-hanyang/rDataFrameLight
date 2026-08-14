#include <stdexcept>
#include <string>
#include <map>
#include <memory>
#include <algorithm>
#include <cmath>
#include <vector>

std::shared_ptr<const correction::CorrectionSet> BTagCorrShape_cset;
std::shared_ptr<const correction::Correction> BTagCorrShape_corr;

std::shared_ptr<const correction::Correction> getBTagShapeCorrectionOrThrow(const std::string& name)
{
    try {
        return BTagCorrShape_cset->at(name);
    } catch (const std::exception& e) {
        throw std::runtime_error("BTag shape correction not found: " + name + " (" + e.what() + ")");
    }
}

std::string get_btag_shape_corr_file(const std::string& era)
{
    std::map<std::string, std::string> corrFile = {
        {"Run3Summer22NanoAODv12", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/BTV/2022_Summer22/btagging.json"},
        {"Run3Summer22EENanoAODv12", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/BTV/2022_Summer22EE/btagging.json"},
        {"Run3Summer23NanoAODv12", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/BTV/2023_Summer23/btagging.json"},
        {"Run3Summer23BPixNanoAODv12", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/BTV/2023_Summer23BPix/btagging.json"},
        {"Run3Summer23BPixNanoAODv12_SSCR", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/BTV/2023_Summer23BPix/btagging.json"},
        {"Run3Summer23BPixNanoAODv12_AR", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/BTV/2023_Summer23BPix/btagging.json"},
    };
    auto fileIt = corrFile.find(era);
    if (fileIt == corrFile.end()) {
        throw std::runtime_error("Unknown era for BTag shape correction file: " + era);
    }
    return fileIt->second;
}

std::string get_btag_shape_payload_name(const std::string& era)
{
    std::map<std::string, std::string> correctionTab{
        {"Run3Summer22NanoAODv12", "particleNet_shape"},
        {"Run3Summer22EENanoAODv12", "particleNet_shape"},
        {"Run3Summer23NanoAODv12", "particleNet_shape"},
        {"Run3Summer23BPixNanoAODv12", "particleNet_shape"},
        {"Run3Summer23BPixNanoAODv12_SSCR", "particleNet_shape"},
        {"Run3Summer23BPixNanoAODv12_AR", "particleNet_shape"},
    };
    auto it = correctionTab.find(era);
    if (it == correctionTab.end()) {
        throw std::runtime_error("Unknown era for btag shape payload: " + era);
    }
    return it->second;
}

void BTagShape_init(const std::string& era)
{
    BTagCorrShape_cset = correction::CorrectionSet::from_file(get_btag_shape_corr_file(era));
    BTagCorrShape_corr = getBTagShapeCorrectionOrThrow(get_btag_shape_payload_name(era));
}

bool btagShapeSystAppliesToFlavor(const std::string& source, unsigned char flavor)
{
    if (source == "cferr1" || source == "cferr2")
        return flavor == 4;
    if (source == "lf" || source == "lfstats1" || source == "lfstats2" ||
        source == "hf" || source == "hfstats1" || source == "hfstats2")
        return flavor == 0 || flavor == 5;
    throw std::runtime_error("Unknown btag shape systematic source: " + source);
}

ROOT::VecOps::RVec<float> get_SF_shape(
    std::string syst,
    const ROOT::VecOps::RVec<unsigned char>& Jet_fla,
    const ROOT::VecOps::RVec<float>& Jet_eta,
    const ROOT::VecOps::RVec<float>& Jet_pt,
    const ROOT::VecOps::RVec<float>& Jet_btagPNetB,
    const ROOT::VecOps::RVec<int>& Jet_passJetIdTightLepVeto)
{
    ROOT::VecOps::RVec<float> Jet_particleNet_shape_nom(Jet_eta.size(), 1.0f);
    for (size_t i = 0; i < Jet_eta.size(); i++)
    {
        if (std::abs(Jet_eta[i]) > 2.5f)
            continue;
        if (!Jet_passJetIdTightLepVeto[i])
            continue;
        if (Jet_btagPNetB[i] < 0.0f)
            continue;

        float absEta = std::abs(std::clamp(Jet_eta[i], -2.49f, 2.49f));
        float safePt = std::clamp(Jet_pt[i], 20.1f, 199.9f);
        float safeScore = std::clamp(Jet_btagPNetB[i], 0.0001f, 1.0999f);
        Jet_particleNet_shape_nom[i] = BTagCorrShape_corr->evaluate(
            std::vector<correction::Variable::Type>{syst, Jet_fla[i], absEta, safePt, safeScore}
        );
    }
    return Jet_particleNet_shape_nom;
}

ROOT::VecOps::RVec<float> get_SF_shape_flavour_syst(
    const std::string& source,
    const std::string& direction,
    const ROOT::VecOps::RVec<unsigned char>& Jet_fla,
    const ROOT::VecOps::RVec<float>& Jet_eta,
    const ROOT::VecOps::RVec<float>& Jet_pt,
    const ROOT::VecOps::RVec<float>& Jet_btagPNetB,
    const ROOT::VecOps::RVec<int>& Jet_passJetIdTightLepVeto)
{
    if (direction != "up" && direction != "down")
        throw std::runtime_error("Unknown btag shape systematic direction: " + direction);

    ROOT::VecOps::RVec<float> Jet_particleNet_shape(Jet_eta.size(), 1.0f);
    for (size_t i = 0; i < Jet_eta.size(); i++)
    {
        if (std::abs(Jet_eta[i]) > 2.5f)
            continue;
        if (!Jet_passJetIdTightLepVeto[i])
            continue;
        if (Jet_btagPNetB[i] < 0.0f)
            continue;

        const std::string syst = btagShapeSystAppliesToFlavor(source, Jet_fla[i]) ? direction + "_" + source : "central";
        float absEta = std::abs(std::clamp(Jet_eta[i], -2.49f, 2.49f));
        float safePt = std::clamp(Jet_pt[i], 20.1f, 199.9f);
        float safeScore = std::clamp(Jet_btagPNetB[i], 0.0001f, 1.0999f);
        Jet_particleNet_shape[i] = BTagCorrShape_corr->evaluate(
            std::vector<correction::Variable::Type>{syst, Jet_fla[i], absEta, safePt, safeScore}
        );
    }
    return Jet_particleNet_shape;
}

float compute_total_weight_shape(const ROOT::VecOps::RVec<float>& Jet_btag_SF, const ROOT::VecOps::RVec<int>& GoodJetCond)
{
    float product = 1.0f;
    for (size_t i = 0; i < Jet_btag_SF.size(); i++)
    {
        if (!GoodJetCond[i]) continue;
        product *= Jet_btag_SF[i];
    }
    return product;
}
