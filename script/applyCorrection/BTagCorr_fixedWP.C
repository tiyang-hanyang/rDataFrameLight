#include <stdexcept>
#include <string>
#include <map>
#include <memory>
#include <algorithm>
#include <cmath>
#include <vector>

std::shared_ptr<const correction::CorrectionSet> BTagCorrFixedWP_cset;
std::shared_ptr<const correction::Correction> BTagCorrFixedWP_corr;
std::shared_ptr<const correction::CorrectionSet> BTagCorrFixedWP_cseteff;
std::shared_ptr<const correction::Correction> BTagCorrFixedWP_ceff;

std::shared_ptr<const correction::Correction> getBTagFixedWPCorrectionOrThrow(const std::string& name)
{
    try {
        return BTagCorrFixedWP_cset->at(name);
    } catch (const std::exception& e) {
        throw std::runtime_error("BTag fixedWP correction not found: " + name + " (" + e.what() + ")");
    }
}

float get_2024_wp_threshold(const std::string& working_point)
{
    if (working_point == "L") return 0.0246f;
    if (working_point == "M") return 0.1272f;
    if (working_point == "T") return 0.4648f;
    throw std::runtime_error("Unknown 2024 working point: " + working_point);
}

std::string get_btag_fixedwp_corr_file(const std::string& era)
{
    std::map<std::string, std::string> corrFile = {
        {"RunIII2024Summer24NanoAODv15", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/BTV/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/btagging.json"},
        {"RunIII2024Summer24NanoAODv15_SSCR", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/BTV/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/btagging.json"},
        {"RunIII2024Summer24NanoAODv15_AR", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/BTV/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/btagging.json"},
    };
    auto fileIt = corrFile.find(era);
    if (fileIt == corrFile.end()) {
        throw std::runtime_error("Unknown era for BTag fixedWP correction file: " + era);
    }
    return fileIt->second;
}

std::string get_btag_fixedwp_payload_name(const std::string& era)
{
    std::map<std::string, std::string> correctionTab{
        {"RunIII2024Summer24NanoAODv15", "UParTAK4_comb"},
        {"RunIII2024Summer24NanoAODv15_SSCR", "UParTAK4_comb"},
        {"RunIII2024Summer24NanoAODv15_AR", "UParTAK4_comb"}
    };
    auto it = correctionTab.find(era);
    if (it == correctionTab.end()) {
        throw std::runtime_error("Unknown era for btag fixedWP payload: " + era);
    }
    return it->second;
}

std::string get_btag_eff_file(const std::string& era, const std::string& channel)
{
    std::string btagEffDir = "/home/tiyang/public/rDataFrameLight_update/source/json/samples/FourJet_NanoAOD/btag_efficiency_"+era+"/";
    if (era == "RunIII2024Summer24NanoAODv15_SSCR")
    {
        btagEffDir = "/home/tiyang/public/rDataFrameLight_update/source/json/samples/SameSign_CR/btag_eff_RunIII2024Summer24NanoAODv15_SSCR/";
    }
    if (era == "RunIII2024Summer24NanoAODv15_AR")
    {
        btagEffDir = "/home/tiyang/public/rDataFrameLight_update/source/json/samples/special_relax_miniisolation/btag_special_relax_miniIso_4jets/";
    }
    return btagEffDir+"/"+channel+"_btag_eff.json";
}

void BtagFixedWP_init(const std::string& era, const std::string& channel)
{
    BTagCorrFixedWP_cset = correction::CorrectionSet::from_file(get_btag_fixedwp_corr_file(era));
    BTagCorrFixedWP_corr = getBTagFixedWPCorrectionOrThrow(get_btag_fixedwp_payload_name(era));

    const std::string effFile = get_btag_eff_file(era, channel);
    BTagCorrFixedWP_cseteff = correction::CorrectionSet::from_file(effFile);
    try {
        BTagCorrFixedWP_ceff = BTagCorrFixedWP_cseteff->at("UParTAK4_eff_values");
    } catch (const std::exception& e) {
        throw std::runtime_error("BTag efficiency correction not found in " + effFile + ": UParTAK4_eff_values (" + e.what() + ")");
    }
}

ROOT::VecOps::RVec<float> get_SF_fixedWP(
    std::string syst,
    std::string working_point,
    const ROOT::VecOps::RVec<unsigned char>& Jet_fla,
    const ROOT::VecOps::RVec<float>& Jet_eta,
    const ROOT::VecOps::RVec<float>& Jet_pt)
{
    ROOT::VecOps::RVec<float> Jet_fixedWP_SF(Jet_eta.size());
    for (size_t i = 0; i < Jet_eta.size(); i++)
    {
        if (Jet_fla[i] == 0) {
            Jet_fixedWP_SF[i] = 1.0f;
            continue;
        }

        float absEta = std::abs(std::clamp(Jet_eta[i], -2.49f, 2.49f));
        float safePt = std::clamp(Jet_pt[i], 20.1f, 999.9f);
        Jet_fixedWP_SF[i] = BTagCorrFixedWP_corr->evaluate(
            std::vector<correction::Variable::Type>{syst, working_point, Jet_fla[i], absEta, safePt}
        );
    }
    return Jet_fixedWP_SF;
}

float get_fixedWP_SF_2024_or_unity_for_light(
    const std::string& syst,
    const std::string& working_point,
    unsigned char flavor,
    float eta,
    float pt)
{
    if (flavor == 0)
        return 1.0f;

    float absEta = std::abs(std::clamp(eta, -2.49f, 2.49f));
    float safePt = std::clamp(pt, 20.1f, 999.9f);
    return BTagCorrFixedWP_corr->evaluate(std::vector<correction::Variable::Type>{syst, working_point, flavor, absEta, safePt});
}

float get_eff_value(
    const std::string& syst,
    const std::string& working_point,
    unsigned char flavor,
    float eta,
    float pt)
{
    float absEta = std::abs(std::clamp(eta, -2.49f, 2.49f));
    float safePt = std::clamp(pt, 20.1f, 999.9f);
    return BTagCorrFixedWP_ceff->evaluate(std::vector<correction::Variable::Type>{syst, working_point, flavor, absEta, safePt});
}

ROOT::VecOps::RVec<float> get_eff(
    std::string syst,
    std::string working_point,
    const ROOT::VecOps::RVec<unsigned char>& Jet_fla,
    const ROOT::VecOps::RVec<float>& Jet_eta,
    const ROOT::VecOps::RVec<float>& Jet_pt)
{
    ROOT::VecOps::RVec<float> Jet_btag_eff(Jet_eta.size());
    for (size_t i = 0; i < Jet_eta.size(); i++)
    {
        float absEta = std::abs(std::clamp(Jet_eta[i], -2.49f, 2.49f));
        float safePt = std::clamp(Jet_pt[i], 20.1f, 999.9f);
        Jet_btag_eff[i] = BTagCorrFixedWP_ceff->evaluate(
            std::vector<correction::Variable::Type>{syst, working_point, Jet_fla[i], absEta, safePt}
        );
    }
    return Jet_btag_eff;
}

float compute_total_weight_fixedWP_2024(
    const ROOT::VecOps::RVec<float>& Jet_btag_SF,
    const ROOT::VecOps::RVec<float>& Jet_btag_eff,
    const ROOT::VecOps::RVec<float>& Jet_btag_score,
    const ROOT::VecOps::RVec<int>& GoodJetCond,
    float threshold)
{
    float product = 1.0f;
    for (size_t i = 0; i < Jet_btag_SF.size(); i++)
    {
        if (!GoodJetCond[i]) continue;

        if (Jet_btag_score[i] > threshold)
        {
            product *= Jet_btag_SF[i];
        }
        else if (Jet_btag_eff[i] < 1.0f)
        {
            product *= (1.0f - Jet_btag_SF[i] * Jet_btag_eff[i]) / (1.0f - Jet_btag_eff[i]);
        }
    }
    return product;
}

float compute_total_weight_fixedWP_2024_syst(
    const std::string& syst,
    const std::string& working_point,
    const ROOT::VecOps::RVec<unsigned char>& Jet_fla,
    const ROOT::VecOps::RVec<float>& Jet_eta,
    const ROOT::VecOps::RVec<float>& Jet_pt,
    const ROOT::VecOps::RVec<float>& Jet_btag_score,
    const ROOT::VecOps::RVec<int>& GoodJetCond,
    float threshold)
{
    float product = 1.0f;
    for (size_t i = 0; i < Jet_fla.size(); i++)
    {
        if (!GoodJetCond[i])
            continue;

        const float sf = get_fixedWP_SF_2024_or_unity_for_light(syst, working_point, Jet_fla[i], Jet_eta[i], Jet_pt[i]);
        if (sf == 1.0f)
            continue;

        const float eff = get_eff_value("central", working_point, Jet_fla[i], Jet_eta[i], Jet_pt[i]);
        if (Jet_btag_score[i] > threshold)
        {
            product *= sf;
        }
        else if (eff < 1.0f)
        {
            product *= (1.0f - sf * eff) / (1.0f - eff);
        }
    }
    return product;
}

float btagWeightQuadratureUp(float nominal, const std::vector<float>& variations)
{
    float sum2 = 0.0f;
    for (const auto& varied : variations)
    {
        const float diff = varied - nominal;
        if (diff > 0.0f)
            sum2 += diff * diff;
    }
    return nominal + std::sqrt(sum2);
}

float btagWeightQuadratureDown(float nominal, const std::vector<float>& variations)
{
    float sum2 = 0.0f;
    for (const auto& varied : variations)
    {
        const float diff = nominal - varied;
        if (diff > 0.0f)
            sum2 += diff * diff;
    }
    return std::max(0.0f, nominal - std::sqrt(sum2));
}
