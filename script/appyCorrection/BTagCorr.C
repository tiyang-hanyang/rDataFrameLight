#include <functional>
#include <stdexcept>
#include <string>
#include <map>
#include <memory>
#include <algorithm>
#include <cmath>

std::shared_ptr<const correction::CorrectionSet> BTagCorr_cset;
std::shared_ptr<const correction::Correction> BTagCorr_corr;
std::shared_ptr<const correction::CorrectionSet> BTagCorr_cseteff;
std::shared_ptr<const correction::Correction> BTagCorr_ceff;

void Btag_init(const std::string& era)
{
    // hard coding the path of the file here (as just macro)
    std::map<std::string, std::string> corrFile = {
        {"RunIII2024Summer24NanoAODv15", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/BTV/2024_Summer24/btagging_preliminary.json"},
        {"RunIII2024Summer24NanoAODv15_SSCR", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/BTV/2024_Summer24/btagging_preliminary.json"},
        {"Run3Summer23NanoAODv12", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/BTV/2023_Summer23/btagging.json"},
        {"Run3Summer23BPixNanoAODv12", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/BTV/2023_Summer23BPix/btagging.json"},
        {"Run3Summer22NanoAODv12", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/BTV/2022_Summer22/btagging.json"},
        {"Run3Summer22EENanoAODv12", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/BTV/2022_Summer22EE/btagging.json"},
    };
    BTagCorr_cset = correction::CorrectionSet::from_file(corrFile[era]);

    // hard coding the name of the table to load
    std::map<std::string, std::string> correctionTab{
        {"Run3Summer22NanoAODv12", "particleNet_shape"},
        {"Run3Summer22EENanoAODv12", "particleNet_shape"},
        {"Run3Summer23NanoAODv12", "particleNet_shape"},
        {"Run3Summer23BPixNanoAODv12", "particleNet_shape"},
        {"RunIII2024Summer24NanoAODv15", "UParTAK4_kinfit"},
        {"RunIII2024Summer24NanoAODv15_SSCR", "UParTAK4_kinfit"}
    };
    auto it = correctionTab.find(era);
    if (it == correctionTab.end()) {
        throw std::runtime_error("Unknown era: " + era);
    }
    BTagCorr_corr = BTagCorr_cset->at(it->second);
}

void reload_eff(const std::string& era, const std::string& channel)
{
    // change to the channel final for loading the corresponding efficiency factors
    std::string btagEffDir = "/home/tiyang/public/rDataFrameLight_update/source/json/samples/FourJet_NanoAOD/btag_efficiency_"+era+"/";
    if (era == "RunIII2024Summer24NanoAODv15_SSCR")
    {
        btagEffDir = "/home/tiyang/public/rDataFrameLight_update/source/json/samples/SameSign_CR/btag_eff_RunIII2024Summer24NanoAODv15_SSCR/";
    }
    BTagCorr_cseteff = correction::CorrectionSet::from_file(btagEffDir+"/"+channel+"_btag_eff.json");
    BTagCorr_ceff =  BTagCorr_cseteff->at("UParTAK4_eff_values");
}

// syst: "central" for nominal input
// flavour: b=5, c=4, udsg=0
ROOT::VecOps::RVec<float> get_SF_shape(std::string syst, const ROOT::VecOps::RVec<unsigned char>& Jet_fla, const ROOT::VecOps::RVec<float>& Jet_eta, const ROOT::VecOps::RVec<float>& Jet_pt, const ROOT::VecOps::RVec<float>& Jet_btagPNetB)
{
    ROOT::VecOps::RVec<float> Jet_particleNet_shape_nom(Jet_eta.size());
    for (size_t i = 0; i < Jet_eta.size(); i++)
    {
        // eta should already in range -2.5-2.5 based on jetvetomap
        float absEta = std::abs(std::clamp(Jet_eta[i], -2.49f, 2.49f));
        float safePt = std::clamp(Jet_pt[i], 20.1f, 199.9f);
        float safeScore = std::clamp(Jet_btagPNetB[i], 0.0001f, 1.0999f);
        Jet_particleNet_shape_nom[i] = BTagCorr_corr->evaluate(std::vector<correction::Variable::Type>{syst, Jet_fla[i], absEta, safePt, safeScore});
    }
    return Jet_particleNet_shape_nom;
}

ROOT::VecOps::RVec<float> get_SF_fixedWP(std::string syst, const ROOT::VecOps::RVec<unsigned char>& Jet_fla, const ROOT::VecOps::RVec<float>& Jet_eta, const ROOT::VecOps::RVec<float>& Jet_pt)
{
    ROOT::VecOps::RVec<float> Jet_fixedWP_SF(Jet_eta.size());
    for (size_t i = 0; i < Jet_eta.size(); i++)
    {
        // eta should already in range -2.5-2.5 based on jetvetomap
        float absEta = std::abs(std::clamp(Jet_eta[i], -2.49f, 2.49f));
        float safePt = std::clamp(Jet_pt[i], 20.1f, 999.9f);
        Jet_fixedWP_SF[i] = BTagCorr_corr->evaluate(std::vector<correction::Variable::Type>{syst, "M", 5, absEta, safePt});
    }
    return Jet_fixedWP_SF;
}

// need the efficiency loading
// I still preseve the API for the systematics, and possible eta splitting
// Now I only need the pt partition and the truth flavour, syst only has central now, absEta only need to be in range.
ROOT::VecOps::RVec<float> get_eff(std::string syst, const ROOT::VecOps::RVec<unsigned char>& Jet_fla, const ROOT::VecOps::RVec<float>& Jet_eta, const ROOT::VecOps::RVec<float>& Jet_pt)
{
    ROOT::VecOps::RVec<float> Jet_btag_eff(Jet_eta.size());
    for (size_t i = 0; i < Jet_eta.size(); i++)
    {
        // eta should already in range -2.5-2.5 based on jetvetomap
        float absEta = std::abs(std::clamp(Jet_eta[i], -2.49f, 2.49f));
        float safePt = std::clamp(Jet_pt[i], 20.1f, 999.9f);
        Jet_btag_eff[i] = BTagCorr_ceff->evaluate(std::vector<correction::Variable::Type>{syst, "M", Jet_fla[i], absEta, safePt});
    }
    return Jet_btag_eff;
}

float compute_total_weight_old(const ROOT::VecOps::RVec<float>& Jet_btag_weight, const ROOT::VecOps::RVec<int>& pass_btag)
{
    float product=1.0;
    for (size_t i = 0; i < pass_btag.size(); i++)
    {
        if (pass_btag[i]==1) product *= Jet_btag_weight[i];
    }
    return product;
}

float compute_total_weight_2024(const ROOT::VecOps::RVec<float>& Jet_btag_SF, const ROOT::VecOps::RVec<float>& Jet_btag_eff, const ROOT::VecOps::RVec<float>& Jet_btag_score, const ROOT::VecOps::RVec<int>& GoodJetCond)
{
    float product=1.0;
    for (size_t i = 0; i < Jet_btag_SF.size(); i++)
    {
        // only need to correct the good jets I used in this study
        if (! GoodJetCond[i]) continue;

        // now hard-coding the medium working points for the Jet_btagUParTAK4B
        if (Jet_btag_score[i] > 0.1272)
        {
            // if pass the btag, then just need the SF = SF * eff / eff
            product *= Jet_btag_SF[i];
        }
        else 
        {
            // if not pass the btag, then (1 - SF * eff) / （1 - eff）
            if (Jet_btag_eff[i]<1)
                product *= (1 -  Jet_btag_SF[i] * Jet_btag_eff[i]) / (1 - Jet_btag_eff[i]);
        }
    }
    return product;
}


float compute_total_weight_2023(const ROOT::VecOps::RVec<float>& Jet_btag_SF, const ROOT::VecOps::RVec<float>& Jet_btag_eff, const ROOT::VecOps::RVec<float>& Jet_btag_score, const ROOT::VecOps::RVec<int>& GoodJetCond)
{
    float product=1.0;
    for (size_t i = 0; i < Jet_btag_SF.size(); i++)
    {
        // only need to correct the good jets I used in this study
        if (! GoodJetCond[i]) continue;

        // now hard-coding the medium working points for the Jet_btagUParTAK4B
        if (Jet_btag_score[i] > 0.1917)
        {
            // if pass the btag, then just need the SF = SF * eff / eff
            product *= Jet_btag_SF[i];
        }
        else 
        {
            // if not pass the btag, then (1 - SF * eff) / （1 - eff）
            if (Jet_btag_eff[i]<1)
                product *= (1 -  Jet_btag_SF[i] * Jet_btag_eff[i]) / (1 - Jet_btag_eff[i]);
        }
    }
    return product;
}


float compute_total_weight_2023BPix(const ROOT::VecOps::RVec<float>& Jet_btag_SF, const ROOT::VecOps::RVec<float>& Jet_btag_eff, const ROOT::VecOps::RVec<float>& Jet_btag_score, const ROOT::VecOps::RVec<int>& GoodJetCond)
{
    float product=1.0;
    for (size_t i = 0; i < Jet_btag_SF.size(); i++)
    {
        // only need to correct the good jets I used in this study
        if (! GoodJetCond[i]) continue;

        // now hard-coding the medium working points for the Jet_btagUParTAK4B
        if (Jet_btag_score[i] > 0.1919)
        {
            // if pass the btag, then just need the SF = SF * eff / eff
            product *= Jet_btag_SF[i];
        }
        else 
        {
            // if not pass the btag, then (1 - SF * eff) / （1 - eff）
            if (Jet_btag_eff[i]<1)
                product *= (1 -  Jet_btag_SF[i] * Jet_btag_eff[i]) / (1 - Jet_btag_eff[i]);
        }
    }
    return product;
}



float compute_total_weight_2022(const ROOT::VecOps::RVec<float>& Jet_btag_SF, const ROOT::VecOps::RVec<float>& Jet_btag_eff, const ROOT::VecOps::RVec<float>& Jet_btag_score, const ROOT::VecOps::RVec<int>& GoodJetCond)
{
    float product=1.0;
    for (size_t i = 0; i < Jet_btag_SF.size(); i++)
    {
        // only need to correct the good jets I used in this study
        if (! GoodJetCond[i]) continue;

        // now hard-coding the medium working points for the Jet_btagUParTAK4B
        if (Jet_btag_score[i] > 0.245)
        {
            // if pass the btag, then just need the SF = SF * eff / eff
            product *= Jet_btag_SF[i];
        }
        else 
        {
            // if not pass the btag, then (1 - SF * eff) / （1 - eff）
            if (Jet_btag_eff[i]<1)
                product *= (1 -  Jet_btag_SF[i] * Jet_btag_eff[i]) / (1 - Jet_btag_eff[i]);
        }
    }
    return product;
}


float compute_total_weight_2022EE(const ROOT::VecOps::RVec<float>& Jet_btag_SF, const ROOT::VecOps::RVec<float>& Jet_btag_eff, const ROOT::VecOps::RVec<float>& Jet_btag_score, const ROOT::VecOps::RVec<int>& GoodJetCond)
{
    float product=1.0;
    for (size_t i = 0; i < Jet_btag_SF.size(); i++)
    {
        // only need to correct the good jets I used in this study
        if (! GoodJetCond[i]) continue;

        // now hard-coding the medium working points for the Jet_btagUParTAK4B
        if (Jet_btag_score[i] > 0.2605)
        {
            // if pass the btag, then just need the SF = SF * eff / eff
            product *= Jet_btag_SF[i];
        }
        else 
        {
            // if not pass the btag, then (1 - SF * eff) / （1 - eff）
            if (Jet_btag_eff[i]<1)
                product *= (1 -  Jet_btag_SF[i] * Jet_btag_eff[i]) / (1 - Jet_btag_eff[i]);
        }
    }
    return product;
}