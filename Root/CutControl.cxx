#include "CutControl.h"
#include "Utility.h"

#include <limits>
#include <stdexcept>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <set>
#include <initializer_list>
#include <cmath>
#include <regex>

#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"
#include "TInterpreter.h"

//////////////////////////////////////////////////
// Free helper functions for the cling calling
//////////////////////////////////////////////////

int getMatchIdx(float eta, float phi,
                const ROOT::VecOps::RVec<float> &etas,
                const ROOT::VecOps::RVec<float> &phis)
{
    int   bestIdx = -1;
    float bestDR  = 0.4f;
    const auto n = etas.size();
    for (size_t i = 0; i < n; ++i) {
        float dr = ROOT::VecOps::DeltaR(eta, etas[i], phi, phis[i]);
        if (dr < bestDR) {
            bestDR  = dr;
            bestIdx = static_cast<int>(i);
        }
    }
    return bestIdx;
}

bool hasOSLooseMuonPairOnZ(
    const ROOT::VecOps::RVec<int> &looseIdx,
    const ROOT::VecOps::RVec<float> &Muon_pt,
    const ROOT::VecOps::RVec<float> &Muon_eta,
    const ROOT::VecOps::RVec<float> &Muon_phi,
    const ROOT::VecOps::RVec<float> &Muon_mass,
    const ROOT::VecOps::RVec<int> &Muon_charge,
    float zMass = 91.2f,
    float window = 10.0f)
{
    for (size_t i = 0; i < looseIdx.size(); ++i)
    {
        for (size_t j = i + 1; j < looseIdx.size(); ++j)
        {
            const int idx1 = looseIdx[i];
            const int idx2 = looseIdx[j];

            if (idx1 < 0 || idx2 < 0)
                continue;
            if (Muon_charge[idx1] * Muon_charge[idx2] >= 0)
                continue;

            TLorentzVector p1, p2;
            p1.SetPtEtaPhiM(Muon_pt[idx1], Muon_eta[idx1], Muon_phi[idx1], Muon_mass[idx1]);
            p2.SetPtEtaPhiM(Muon_pt[idx2], Muon_eta[idx2], Muon_phi[idx2], Muon_mass[idx2]);

            if (std::abs((p1 + p2).M() - zMass) < window)
                return true;
        }
    }
    return false;
}

namespace {
    std::string escapeRegex(const std::string &text)
    {
        static const std::regex re{R"([-[\]{}()*+?.,\^$|#\s])"};
        return std::regex_replace(text, re, R"(\$&)");
    }

    std::string replaceIdentifierTokens(
        const std::string &expression,
        const std::map<std::string, std::string> &replacements)
    {
        std::string updated = expression;
        std::vector<std::pair<std::string, std::string>> ordered(replacements.begin(), replacements.end());
        std::sort(
            ordered.begin(),
            ordered.end(),
            [](const auto &a, const auto &b) { return a.first.size() > b.first.size(); });

        for (const auto &[from, to] : ordered)
        {
            const std::regex tokenRegex("(^|[^A-Za-z0-9_])(" + escapeRegex(from) + ")(?=[^A-Za-z0-9_]|$)");
            updated = std::regex_replace(updated, tokenRegex, "$1" + to);
        }
        return updated;
    }

    std::string sanitizeCutTagName(const std::string &name)
    {
        std::string sanitized;
        sanitized.reserve(name.size());
        bool lastUnderscore = false;
        for (char ch : name)
        {
            const bool keep = (ch >= 'A' && ch <= 'Z')
                           || (ch >= 'a' && ch <= 'z')
                           || (ch >= '0' && ch <= '9');
            if (keep)
            {
                sanitized.push_back(ch);
                lastUnderscore = false;
            }
            else if (!lastUnderscore)
            {
                sanitized.push_back('_');
                lastUnderscore = true;
            }
        }
        while (!sanitized.empty() && sanitized.back() == '_')
            sanitized.pop_back();
        if (sanitized.empty())
            return "unnamed_cut";
        return sanitized;
    }

    void DeclareHelpersToCling()
    {
        static bool done = false;
        if (done) return;
        gInterpreter->Declare(R"(
            #include "ROOT/RVec.hxx"
            #include <algorithm>
            #include <cmath>
            #include <initializer_list>
            #include <stdexcept>
            int getMatchIdx(float eta, float phi,
                const ROOT::VecOps::RVec<float> &etas,
                const ROOT::VecOps::RVec<float> &phis)
            {
                int   bestIdx = -1;
                float bestDR  = 0.4f;
                const auto n = etas.size();
                for (size_t i = 0; i < n; ++i) {
                    float dr = ROOT::VecOps::DeltaR(eta, etas[i], phi, phis[i]);
                    if (dr < bestDR) {
                        bestDR  = dr;
                        bestIdx = static_cast<int>(i);
                    }
                }
                return bestIdx;
            }
            bool hasOSLooseMuonPairOnZ(
                const ROOT::VecOps::RVec<int> &looseIdx,
                const ROOT::VecOps::RVec<float> &Muon_pt,
                const ROOT::VecOps::RVec<float> &Muon_eta,
                const ROOT::VecOps::RVec<float> &Muon_phi,
                const ROOT::VecOps::RVec<float> &Muon_mass,
                const ROOT::VecOps::RVec<int> &Muon_charge,
                float zMass = 91.2f,
                float window = 10.0f)
            {
                for (size_t i = 0; i < looseIdx.size(); ++i)
                {
                    for (size_t j = i + 1; j < looseIdx.size(); ++j)
                    {
                        const int idx1 = looseIdx[i];
                        const int idx2 = looseIdx[j];

                        if (idx1 < 0 || idx2 < 0)
                            continue;
                        if (Muon_charge[idx1] * Muon_charge[idx2] >= 0)
                            continue;

                        TLorentzVector p1, p2;
                        p1.SetPtEtaPhiM(Muon_pt[idx1], Muon_eta[idx1], Muon_phi[idx1], Muon_mass[idx1]);
                        p2.SetPtEtaPhiM(Muon_pt[idx2], Muon_eta[idx2], Muon_phi[idx2], Muon_mass[idx2]);

                        if (std::abs((p1 + p2).M() - zMass) < window)
                            return true;
                    }
                }
                return false;
            }
            ROOT::VecOps::RVec<float> allJetDRFromMuon( ROOT::VecOps::RVec<float>Jet_eta, float Muon_eta,  ROOT::VecOps::RVec<float> Jet_phi, float Muon_phi) {
                 ROOT::VecOps::RVec<float> DR;
                int nJet = Jet_eta.size();
                for (int i =0;i <nJet; i++) {
                    DR.push_back(ROOT::VecOps::DeltaR(Jet_eta[i], Muon_eta, Jet_phi[i], Muon_phi));
                }
                return DR;
            }
            ROOT::VecOps::RVec<float> getMuonMatchJetDR(const ROOT::VecOps::RVec<short>& muon_jetIdx, const ROOT::VecOps::RVec<float>& jet_eta, const ROOT::VecOps::RVec<float>& jet_phi, const ROOT::VecOps::RVec<float>& muon_eta, const ROOT::VecOps::RVec<float>& muon_phi) 
            {
                ROOT::VecOps::RVec<float> muonjetdr;
                int nMuon = muon_jetIdx.size();
                for (int i=0; i < nMuon; i++) {
                    if(muon_jetIdx[i]==-1) {
                        muonjetdr.push_back(999);
                        continue;
                    }
                    muonjetdr.push_back(ROOT::VecOps::DeltaR(muon_eta[i], jet_eta[muon_jetIdx[i]], muon_phi[i], jet_phi[muon_jetIdx[i]]));
                }
                return muonjetdr;
            }
            ROOT::VecOps::RVec<float> getMuonMatchJetPt(const ROOT::VecOps::RVec<short>& muon_jetIdx, const ROOT::VecOps::RVec<float>& jet_pt) 
            {
                ROOT::VecOps::RVec<float> muonjetpt;
                int nMuon = muon_jetIdx.size();
                for (int i=0; i < nMuon; i++) {
                    if(muon_jetIdx[i]==-1) {
                        muonjetpt.push_back(-999);
                        continue;
                    }
                    muonjetpt.push_back(jet_pt[muon_jetIdx[i]]);
                }
                return muonjetpt;
            }
            ROOT::VecOps::RVec<short> ensureMuonJetGood(const ROOT::VecOps::RVec<short>& muon_jetIdx, const ROOT::VecOps::RVec<bool>& GoodJetCond) 
            {
                ROOT::VecOps::RVec<short> Muon_jetIdxGood;
                int nMuon = muon_jetIdx.size();
                for (int i=0; i < nMuon; i++) {
                    if(muon_jetIdx[i]==-1) {
                        Muon_jetIdxGood.push_back(-1);
                        continue;
                    }
                    if (GoodJetCond[muon_jetIdx[i]])
                    {
                        Muon_jetIdxGood.push_back(muon_jetIdx[i]);
                    }
                    else Muon_jetIdxGood.push_back(-1);
                }
                return Muon_jetIdxGood;
            }
            int hadronFlavorFromPdg(int pdgId)
            {
                int apid = pdgId < 0 ? -pdgId : pdgId;
                // heavy-flavor hadrons: 4xx/5xx (mesons), 4xxx/5xxx (baryons), 4xxxx/5xxxx (excited)
                if ((apid >= 500 && apid < 600) || (apid >= 5000 && apid < 6000) || (apid >= 50000 && apid < 60000))
                    return 5;
                if ((apid >= 400 && apid < 500) || (apid >= 4000 && apid < 5000) || (apid >= 40000 && apid < 50000))
                    return 4;
                return 0;
            }
            int classifyOriginFromGen(int gpIdx,
                                      const ROOT::VecOps::RVec<int>& GenPart_pdgId,
                                      const ROOT::VecOps::RVec<int>& GenPart_genPartIdxMother)
            {
                if (gpIdx < 0) return 0;
                const int n = GenPart_pdgId.size();
                int idx = gpIdx;
                int maxDepth = 100;
                int steps = 0;
                int hasC = 0;
                while (idx >= 0 && idx < n && steps < maxDepth) {
                    int flav = hadronFlavorFromPdg(GenPart_pdgId[idx]);
                    if (flav == 5) return 5;
                    if (flav == 4) hasC = 1;
                    idx = GenPart_genPartIdxMother[idx];
                    steps++;
                }
                if (hasC) return 4;
                return 3;
            }
            int classifyOriginFromGenWithStatus(int gpIdx,
                                                const ROOT::VecOps::RVec<int>& GenPart_pdgId,
                                                const ROOT::VecOps::RVec<int>& GenPart_genPartIdxMother,
                                                const ROOT::VecOps::RVec<int>& GenPart_statusFlags)
            {
                if (gpIdx < 0) return 0;
                const int n = GenPart_pdgId.size();
                if (gpIdx >= n) return 0;

                int status = GenPart_statusFlags[gpIdx];
                int isPrompt = (status >> 0) & 1;
                int isPromptTau = ((status >> 5) & 1) || ((status >> 3) & 1);
                if (isPromptTau) return 15;
                if (isPrompt) return 1;

                int idx = gpIdx;
                int maxDepth = 100;
                int steps = 0;
                int hasC = 0;
                while (idx >= 0 && idx < n && steps < maxDepth) {
                    int flav = hadronFlavorFromPdg(GenPart_pdgId[idx]);
                    if (flav == 5) return 5;
                    if (flav == 4) hasC = 1;
                    idx = GenPart_genPartIdxMother[idx];
                    steps++;
                }
                if (hasC) return 4;
                return 3;
            }
            ROOT::VecOps::RVec<int> getOriginFlagsFromGenIdx(
                const ROOT::VecOps::RVec<int>& lepGenIdx,
                const ROOT::VecOps::RVec<int>& GenPart_pdgId,
                const ROOT::VecOps::RVec<int>& GenPart_genPartIdxMother)
            {
                ROOT::VecOps::RVec<int> out;
                out.reserve(lepGenIdx.size());
                for (auto idx : lepGenIdx) {
                    out.push_back(classifyOriginFromGen(idx, GenPart_pdgId, GenPart_genPartIdxMother));
                }
                return out;
            }
            ROOT::VecOps::RVec<int> getOriginFlagsFromGenIdx(
                const ROOT::VecOps::RVec<int>& lepGenIdx,
                const ROOT::VecOps::RVec<int>& GenPart_pdgId,
                const ROOT::VecOps::RVec<int>& GenPart_genPartIdxMother,
                const ROOT::VecOps::RVec<int>& GenPart_statusFlags)
            {
                ROOT::VecOps::RVec<int> out;
                out.reserve(lepGenIdx.size());
                for (auto idx : lepGenIdx) {
                    out.push_back(classifyOriginFromGenWithStatus(idx, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_statusFlags));
                }
                return out;
            }
        )");
        done = true;
    }

    ROOT::RDF::RNode defineOrRedefineExpr(
        ROOT::RDF::RNode rdf,
        const std::string &name,
        const std::string &expression)
    {
        if (rdf.HasColumn(name))
            return rdf.Redefine(name, expression);
        return rdf.Define(name, expression);
    }

    template <typename F, typename... ColNames>
    ROOT::RDF::RNode defineOrRedefineCallable(
        ROOT::RDF::RNode rdf,
        const std::string &name,
        F &&callable,
        ColNames &&...cols)
    {
        if (rdf.HasColumn(name))
            return rdf.Redefine(name, std::forward<F>(callable), std::forward<ColNames>(cols)...);
        return rdf.Define(name, std::forward<F>(callable), std::forward<ColNames>(cols)...);
    }
}


////////////////////////////////////////////////// Constructors

CutControl::CutControl(const std::string &jsonFileName) : _applyLambda(std::nullopt)
{
    this->loadProcedure(jsonFileName);
}

CutControl::CutControl(nlohmann::json json) : _applyLambda(std::nullopt)
{
    this->loadProcedure(json);
}

CutControl::CutControl(const CutControl &origControl) : _steps(origControl._steps), _applyLambda(std::nullopt) {}

////////////////////////////////////////////////// Loading Infor

void CutControl::loadProcedure(const std::string &jsonFileName)
{
    auto jsonSteps = rdfWS_utility::readJson("CutControl", jsonFileName);
    this->loadProcedure(jsonSteps);
}

void CutControl::loadProcedure(nlohmann::json json)
{
    for (const auto &step : json)
        this->_steps.push_back(step);
}

///////////////////////////////////////////// Step Manipulation
void CutControl::addStep(const cut_string moreStep)
{
    this->_steps.push_back(moreStep);
    this->_applyLambda.reset();
}

CutControl CutControl::extend(const CutControl &addCon) const
{
    CutControl newCC(*this);
    for (const auto &step : addCon._steps)
    {
        newCC.addStep(step);
    }
    return newCC;
}

CutControl operator+(const CutControl &con1, const CutControl &con2)
{
    return con1.extend(con2);
}

///////////////////////////////////////////// Invoking
ROOT::RDF::RNode CutControl::applyCut(ROOT::RDF::RNode origRDF)
{
    DeclareHelpersToCling();

    // in case the cut function is not initialized
    if (!(this->_applyLambda))
    {
        this->_applyLambda = [this](ROOT::RDF::RNode origRDF)
        {
            for (const auto &step : this->_steps)
            {
                std::string operation = std::get<0>(step);
                if (operation == "define")
                {
                    origRDF = defineOrRedefineExpr(origRDF, std::get<1>(step), std::get<2>(step));
                }
                else if (operation == "redefine")
                {
                    origRDF = origRDF.Redefine(std::get<1>(step), std::get<2>(step));
                }
                else if (operation == "cut")
                {
                    origRDF = origRDF.Filter(std::get<2>(step));
                }
                else if (operation == "TLVPtEtaPhiM")
                {
                    // parse the TLorentzVector input
                    auto capturedVar = (this->extractTLVComp)(std::get<2>(step));
                    origRDF = defineOrRedefineCallable(origRDF, std::get<1>(step), [](float pt, float eta, float phi, float m)
                                                       {
                        TLorentzVector p4;
                        p4.SetPtEtaPhiM(pt, eta, phi, m);
                        return p4; }, capturedVar);
                }
                else if (operation == "TLVPtEtaPhiM_corr")
                {
                    // parse the TLorentzVector input
                    auto capturedVar = (this->extractTLVComp)(std::get<2>(step));
                    origRDF = defineOrRedefineCallable(origRDF, std::get<1>(step), [](double pt, float eta, float phi, float m)
                                                       {
                        TLorentzVector p4;
                        p4.SetPtEtaPhiM(pt, eta, phi, m);
                        return p4; }, capturedVar);
                }
                else if (operation == "TLVPtEtaPhiE")
                {
                    // parse the TLorentzVector input
                    auto capturedVar = this->extractTLVComp(std::get<2>(step));
                    origRDF = defineOrRedefineCallable(origRDF, std::get<1>(step), [](float pt, float eta, float phi, float e)
                                                       {
                        TLorentzVector p4;
                        p4.SetPtEtaPhiE(pt, eta, phi, e);
                        return p4; }, capturedVar);
                }
                else if (operation == "defineDR")
                {
                    // parse the TLorentzVector input
                    auto capturedVar = this->extractTLVComp(std::get<2>(step));
                    origRDF = defineOrRedefineCallable(origRDF, std::get<1>(step), [](float eta1, float eta2, float phi1, float phi2)
                                                       {
                        auto deta2 = pow(eta1-eta2, 2.0);
                        auto dphi = abs(phi1-phi2);
                        if (dphi > 3.1416) dphi = 6.2832 - dphi;
                        auto dphi2 = pow(dphi, 2.0);
                        auto dr = pow(deta2+dphi2, 0.5);
                        return dr; }, capturedVar);
                }
                else
                {
                    throw std::runtime_error("[CutControl] Operation type of the step not defined: " + operation + ". Must one of define, cut, TLVPtEtaPhiM and TLVPtEtaPhiE! Please check your json file");
                }
            }
            return origRDF;
        };
    }

    return this->_applyLambda.value()(origRDF);
}

ROOT::RDF::RNode CutControl::applyCutSkippingSteps(ROOT::RDF::RNode origRDF, const std::set<std::string> &skipStepNames)
{
    DeclareHelpersToCling();

    for (const auto &step : this->_steps)
    {
        std::string name = std::get<1>(step);
        if (skipStepNames.find(name) != skipStepNames.end())
            continue;

        std::string operation = std::get<0>(step);
        if (operation == "define")
        {
            origRDF = defineOrRedefineExpr(origRDF, name, std::get<2>(step));
        }
        else if (operation == "redefine")
        {
            origRDF = origRDF.Redefine(name, std::get<2>(step));
        }
        else if (operation == "cut")
        {
            origRDF = origRDF.Filter(std::get<2>(step));
        }
        else if (operation == "TLVPtEtaPhiM")
        {
            auto capturedVar = (this->extractTLVComp)(std::get<2>(step));
            origRDF = defineOrRedefineCallable(origRDF, name, [](float pt, float eta, float phi, float m)
                                               {
                TLorentzVector p4;
                p4.SetPtEtaPhiM(pt, eta, phi, m);
                return p4; }, capturedVar);
        }
        else if (operation == "TLVPtEtaPhiM_corr")
        {
            auto capturedVar = (this->extractTLVComp)(std::get<2>(step));
            origRDF = defineOrRedefineCallable(origRDF, name, [](double pt, float eta, float phi, float m)
                                               {
                TLorentzVector p4;
                p4.SetPtEtaPhiM(pt, eta, phi, m);
                return p4; }, capturedVar);
        }
        else if (operation == "TLVPtEtaPhiE")
        {
            auto capturedVar = this->extractTLVComp(std::get<2>(step));
            origRDF = defineOrRedefineCallable(origRDF, name, [](float pt, float eta, float phi, float e)
                                               {
                TLorentzVector p4;
                p4.SetPtEtaPhiE(pt, eta, phi, e);
                return p4; }, capturedVar);
        }
        else if (operation == "defineDR")
        {
            auto capturedVar = this->extractTLVComp(std::get<2>(step));
            origRDF = defineOrRedefineCallable(origRDF, name, [](float eta1, float eta2, float phi1, float phi2)
                                               {
                auto deta2 = pow(eta1-eta2, 2.0);
                auto dphi = abs(phi1-phi2);
                if (dphi > 3.1416) dphi = 6.2832 - dphi;
                auto dphi2 = pow(dphi, 2.0);
                auto dr = pow(deta2+dphi2, 0.5);
                return dr; }, capturedVar);
        }
        else
        {
            throw std::runtime_error("[CutControl] Operation type of the step not defined: " + operation + ". Must one of define, redefine, cut, TLVPtEtaPhiM and TLVPtEtaPhiE! Please check your json file");
        }
    }

    return origRDF;
}

ROOT::RDF::RNode CutControl::applySuffixedCutSubset(
    ROOT::RDF::RNode origRDF,
    const std::set<std::string> &cutStepNamesToApply,
    const std::string &suffix,
    const std::map<std::string, std::string> &baseReplacements)
{
    DeclareHelpersToCling();

    std::map<std::string, std::string> localNames = baseReplacements;
    std::set<std::string> createdNames;
    for (const auto &step : this->_steps)
    {
        const std::string operation = std::get<0>(step);
        const std::string name = std::get<1>(step);

        if (operation == "cut")
        {
            if (cutStepNamesToApply.find(name) == cutStepNamesToApply.end())
                continue;
            const std::string expr = replaceIdentifierTokens(std::get<2>(step), localNames);
            origRDF = origRDF.Filter(expr);
            continue;
        }

        const std::string suffixedName = name + suffix;
        localNames[name] = suffixedName;

        if (operation == "define" || operation == "redefine")
        {
            const std::string expr = replaceIdentifierTokens(std::get<2>(step), localNames);
            if (createdNames.insert(suffixedName).second)
                origRDF = origRDF.Define(suffixedName, expr);
            else
                origRDF = origRDF.Redefine(suffixedName, expr);
        }
        else if (operation == "TLVPtEtaPhiM")
        {
            auto capturedVar = this->extractTLVComp(std::get<2>(step));
            for (auto &var : capturedVar)
            {
                auto found = localNames.find(var);
                if (found != localNames.end())
                    var = found->second;
            }
            if (createdNames.insert(suffixedName).second)
            {
                origRDF = origRDF.Define(suffixedName, [](float pt, float eta, float phi, float m)
                                         {
                TLorentzVector p4;
                p4.SetPtEtaPhiM(pt, eta, phi, m);
                return p4; }, capturedVar);
            }
            else
            {
                origRDF = origRDF.Redefine(suffixedName, [](float pt, float eta, float phi, float m)
                                           {
                TLorentzVector p4;
                p4.SetPtEtaPhiM(pt, eta, phi, m);
                return p4; }, capturedVar);
            }
        }
        else if (operation == "TLVPtEtaPhiM_corr")
        {
            auto capturedVar = this->extractTLVComp(std::get<2>(step));
            for (auto &var : capturedVar)
            {
                auto found = localNames.find(var);
                if (found != localNames.end())
                    var = found->second;
            }
            if (createdNames.insert(suffixedName).second)
            {
                origRDF = origRDF.Define(suffixedName, [](double pt, float eta, float phi, float m)
                                         {
                TLorentzVector p4;
                p4.SetPtEtaPhiM(pt, eta, phi, m);
                return p4; }, capturedVar);
            }
            else
            {
                origRDF = origRDF.Redefine(suffixedName, [](double pt, float eta, float phi, float m)
                                           {
                TLorentzVector p4;
                p4.SetPtEtaPhiM(pt, eta, phi, m);
                return p4; }, capturedVar);
            }
        }
        else if (operation == "TLVPtEtaPhiE")
        {
            auto capturedVar = this->extractTLVComp(std::get<2>(step));
            for (auto &var : capturedVar)
            {
                auto found = localNames.find(var);
                if (found != localNames.end())
                    var = found->second;
            }
            if (createdNames.insert(suffixedName).second)
            {
                origRDF = origRDF.Define(suffixedName, [](float pt, float eta, float phi, float e)
                                         {
                TLorentzVector p4;
                p4.SetPtEtaPhiE(pt, eta, phi, e);
                return p4; }, capturedVar);
            }
            else
            {
                origRDF = origRDF.Redefine(suffixedName, [](float pt, float eta, float phi, float e)
                                           {
                TLorentzVector p4;
                p4.SetPtEtaPhiE(pt, eta, phi, e);
                return p4; }, capturedVar);
            }
        }
        else if (operation == "defineDR")
        {
            auto capturedVar = this->extractTLVComp(std::get<2>(step));
            for (auto &var : capturedVar)
            {
                auto found = localNames.find(var);
                if (found != localNames.end())
                    var = found->second;
            }
            if (createdNames.insert(suffixedName).second)
            {
                origRDF = origRDF.Define(suffixedName, [](float eta1, float eta2, float phi1, float phi2)
                                         {
                auto deta2 = pow(eta1-eta2, 2.0);
                auto dphi = abs(phi1-phi2);
                if (dphi > 3.1416) dphi = 6.2832 - dphi;
                auto dphi2 = pow(dphi, 2.0);
                auto dr = pow(deta2+dphi2, 0.5);
                return dr; }, capturedVar);
            }
            else
            {
                origRDF = origRDF.Redefine(suffixedName, [](float eta1, float eta2, float phi1, float phi2)
                                           {
                auto deta2 = pow(eta1-eta2, 2.0);
                auto dphi = abs(phi1-phi2);
                if (dphi > 3.1416) dphi = 6.2832 - dphi;
                auto dphi2 = pow(dphi, 2.0);
                auto dr = pow(deta2+dphi2, 0.5);
                return dr; }, capturedVar);
            }
        }
        else
        {
            throw std::runtime_error("[CutControl] Operation type of the step not defined: " + operation + ". Must one of define, redefine, cut, TLVPtEtaPhiM and TLVPtEtaPhiE! Please check your json file");
        }
    }

    return origRDF;
}

ROOT::RDF::RNode CutControl::applySuffixedCutTags(
    ROOT::RDF::RNode origRDF,
    const std::string &suffix,
    const std::map<std::string, std::string> &baseReplacements,
    const std::string &finalPassName)
{
    DeclareHelpersToCling();

    std::map<std::string, std::string> localNames = baseReplacements;
    std::set<std::string> createdNames;
    std::string lastPassName;

    for (const auto &step : this->_steps)
    {
        const std::string operation = std::get<0>(step);
        const std::string name = std::get<1>(step);

        if (operation == "cut")
        {
            const std::string expr = replaceIdentifierTokens(std::get<2>(step), localNames);
            const std::string passName = "pass_" + sanitizeCutTagName(name) + suffix;
            const std::string passExpr = (lastPassName.empty() ? expr : "(" + lastPassName + ") && (" + expr + ")");
            if (createdNames.insert(passName).second)
                origRDF = origRDF.Define(passName, passExpr);
            else
                origRDF = origRDF.Redefine(passName, passExpr);
            lastPassName = passName;
            continue;
        }

        const std::string suffixedName = name + suffix;
        localNames[name] = suffixedName;

        if (operation == "define" || operation == "redefine")
        {
            const std::string expr = replaceIdentifierTokens(std::get<2>(step), localNames);
            if (createdNames.insert(suffixedName).second)
                origRDF = origRDF.Define(suffixedName, expr);
            else
                origRDF = origRDF.Redefine(suffixedName, expr);
        }
        else if (operation == "TLVPtEtaPhiM")
        {
            auto capturedVar = this->extractTLVComp(std::get<2>(step));
            for (auto &var : capturedVar)
            {
                auto found = localNames.find(var);
                if (found != localNames.end())
                    var = found->second;
            }
            if (createdNames.insert(suffixedName).second)
            {
                origRDF = origRDF.Define(suffixedName, [](float pt, float eta, float phi, float m)
                                         {
                    TLorentzVector p4;
                    p4.SetPtEtaPhiM(pt, eta, phi, m);
                    return p4; }, capturedVar);
            }
            else
            {
                origRDF = origRDF.Redefine(suffixedName, [](float pt, float eta, float phi, float m)
                                           {
                    TLorentzVector p4;
                    p4.SetPtEtaPhiM(pt, eta, phi, m);
                    return p4; }, capturedVar);
            }
        }
        else if (operation == "TLVPtEtaPhiM_corr")
        {
            auto capturedVar = this->extractTLVComp(std::get<2>(step));
            for (auto &var : capturedVar)
            {
                auto found = localNames.find(var);
                if (found != localNames.end())
                    var = found->second;
            }
            if (createdNames.insert(suffixedName).second)
            {
                origRDF = origRDF.Define(suffixedName, [](double pt, float eta, float phi, float m)
                                         {
                    TLorentzVector p4;
                    p4.SetPtEtaPhiM(pt, eta, phi, m);
                    return p4; }, capturedVar);
            }
            else
            {
                origRDF = origRDF.Redefine(suffixedName, [](double pt, float eta, float phi, float m)
                                           {
                    TLorentzVector p4;
                    p4.SetPtEtaPhiM(pt, eta, phi, m);
                    return p4; }, capturedVar);
            }
        }
        else if (operation == "TLVPtEtaPhiE")
        {
            auto capturedVar = this->extractTLVComp(std::get<2>(step));
            for (auto &var : capturedVar)
            {
                auto found = localNames.find(var);
                if (found != localNames.end())
                    var = found->second;
            }
            if (createdNames.insert(suffixedName).second)
            {
                origRDF = origRDF.Define(suffixedName, [](float pt, float eta, float phi, float e)
                                         {
                    TLorentzVector p4;
                    p4.SetPtEtaPhiE(pt, eta, phi, e);
                    return p4; }, capturedVar);
            }
            else
            {
                origRDF = origRDF.Redefine(suffixedName, [](float pt, float eta, float phi, float e)
                                           {
                    TLorentzVector p4;
                    p4.SetPtEtaPhiE(pt, eta, phi, e);
                    return p4; }, capturedVar);
            }
        }
        else if (operation == "defineDR")
        {
            auto capturedVar = this->extractTLVComp(std::get<2>(step));
            for (auto &var : capturedVar)
            {
                auto found = localNames.find(var);
                if (found != localNames.end())
                    var = found->second;
            }
            if (createdNames.insert(suffixedName).second)
            {
                origRDF = origRDF.Define(suffixedName, [](float eta1, float eta2, float phi1, float phi2)
                                         {
                    auto deta2 = pow(eta1-eta2, 2.0);
                    auto dphi = abs(phi1-phi2);
                    if (dphi > 3.1416) dphi = 6.2832 - dphi;
                    auto dphi2 = pow(dphi, 2.0);
                    auto dr = pow(deta2+dphi2, 0.5);
                    return dr; }, capturedVar);
            }
            else
            {
                origRDF = origRDF.Redefine(suffixedName, [](float eta1, float eta2, float phi1, float phi2)
                                           {
                    auto deta2 = pow(eta1-eta2, 2.0);
                    auto dphi = abs(phi1-phi2);
                    if (dphi > 3.1416) dphi = 6.2832 - dphi;
                    auto dphi2 = pow(dphi, 2.0);
                    auto dr = pow(deta2+dphi2, 0.5);
                    return dr; }, capturedVar);
            }
        }
        else
        {
            throw std::runtime_error("[CutControl] Operation type of the step not defined: " + operation + ". Must one of define, redefine, cut, TLVPtEtaPhiM and TLVPtEtaPhiE! Please check your json file");
        }
    }

    const std::string finalName = finalPassName.empty() ? "pass_all" + suffix : finalPassName;
    if (createdNames.insert(finalName).second)
        origRDF = origRDF.Define(finalName, lastPassName.empty() ? "true" : lastPassName);
    else
        origRDF = origRDF.Redefine(finalName, lastPassName.empty() ? "true" : lastPassName);

    return origRDF;
}

ROOT::RDF::RNode CutControl::applyDefineOnly(
    ROOT::RDF::RNode origRDF,
    const std::function<bool(const std::string&)> &shouldDefine,
    const std::function<void(const std::string&)> &onDefined)
{
    DeclareHelpersToCling();

    for (const auto &step : this->_steps)
    {
        std::string operation = std::get<0>(step);
        if (operation == "cut")
            continue;

        std::string name = std::get<1>(step);
        if (shouldDefine && !shouldDefine(name))
            continue;

        if (operation == "define")
        {
            origRDF = defineOrRedefineExpr(origRDF, name, std::get<2>(step));
        }
        else if (operation == "redefine")
        {
            origRDF = origRDF.Redefine(name, std::get<2>(step));
        }
        else if (operation == "TLVPtEtaPhiM")
        {
            auto capturedVar = (this->extractTLVComp)(std::get<2>(step));
            origRDF = defineOrRedefineCallable(origRDF, name, [](float pt, float eta, float phi, float m)
                                               {
                TLorentzVector p4;
                p4.SetPtEtaPhiM(pt, eta, phi, m);
                return p4; }, capturedVar);
        }
        else if (operation == "TLVPtEtaPhiM_corr")
        {
            auto capturedVar = (this->extractTLVComp)(std::get<2>(step));
            origRDF = defineOrRedefineCallable(origRDF, name, [](double pt, float eta, float phi, float m)
                                               {
                TLorentzVector p4;
                p4.SetPtEtaPhiM(pt, eta, phi, m);
                return p4; }, capturedVar);
        }
        else if (operation == "TLVPtEtaPhiE")
        {
            auto capturedVar = this->extractTLVComp(std::get<2>(step));
            origRDF = defineOrRedefineCallable(origRDF, name, [](float pt, float eta, float phi, float e)
                                               {
                TLorentzVector p4;
                p4.SetPtEtaPhiE(pt, eta, phi, e);
                return p4; }, capturedVar);
        }
        else if (operation == "defineDR")
        {
            auto capturedVar = this->extractTLVComp(std::get<2>(step));
            origRDF = defineOrRedefineCallable(origRDF, name, [](float eta1, float eta2, float phi1, float phi2)
                                               {
                auto deta2 = pow(eta1-eta2, 2.0);
                auto dphi = abs(phi1-phi2);
                if (dphi > 3.1416) dphi = 6.2832 - dphi;
                auto dphi2 = pow(dphi, 2.0);
                auto dr = pow(deta2+dphi2, 0.5);
                return dr; }, capturedVar);
        }
        else
        {
            throw std::runtime_error("[CutControl] Operation type of the step not defined: " + operation + ". Must one of define, cut, TLVPtEtaPhiM and TLVPtEtaPhiE! Please check your json file");
        }

        if (onDefined)
            onDefined(name);
    }

    return origRDF;
}

std::string CutControl::resolveSuffixedColumnName(
    const std::string &baseName,
    const std::string &suffix,
    const std::map<std::string, std::string> &baseReplacements) const
{
    auto replacementIt = baseReplacements.find(baseName);
    if (replacementIt != baseReplacements.end())
        return replacementIt->second;

    for (const auto &step : this->_steps)
    {
        if (std::get<0>(step) == "cut")
            continue;
        if (std::get<1>(step) == baseName)
            return baseName + suffix;
    }

    return baseName;
}

std::vector<std::string> CutControl::extractTLVComp(const std::string &TLVComp) const
{
    std::vector<std::string> totalVars;
    std::string remain = TLVComp;
    while (remain.find(",") != std::string::npos)
    {
        int pos = remain.find(",");
        totalVars.push_back(remain.substr(0, pos));
        remain = remain.substr(pos + 1);
    }
    totalVars.push_back(remain);

    // check if the var is as desired
    if (totalVars.size() != 4)
    {
        throw std::runtime_error("[CutControl] TLorentzVector capture must by 4 variables! Please check your json file and split each by \",\" (no space). Your input is: " + TLVComp);
    }

    return totalVars;
}

///////////////////////////////////////////// Debugging
void CutControl::printSteps()
{
    std::cout << "Printing all the steps" << std::endl;
    for (const auto &step : this->_steps)
    {
        std::cout << "Operation type: " << std::get<0>(step) << "\n\t name: " << std::get<1>(step) << "\n\t expression: " << std::get<2>(step) << std::endl;
    }
}
