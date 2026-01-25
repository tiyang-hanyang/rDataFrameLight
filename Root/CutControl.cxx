#include "CutControl.h"
#include "Utility.h"

#include <limits>
#include <stdexcept>
#include <iostream>
#include <filesystem>
#include <fstream>

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

namespace {
    void DeclareHelpersToCling()
    {
        static bool done = false;
        if (done) return;
        gInterpreter->Declare(R"(
            #include "ROOT/RVec.hxx"
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
        )");
        done = true;
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
                    origRDF = origRDF.Define(std::get<1>(step), std::get<2>(step));
                }
                else if (operation == "cut")
                {
                    origRDF = origRDF.Filter(std::get<2>(step));
                }
                else if (operation == "TLVPtEtaPhiM")
                {
                    // parse the TLorentzVector input
                    auto capturedVar = (this->extractTLVComp)(std::get<2>(step));
                    origRDF = origRDF.Define(std::get<1>(step), [](float pt, float eta, float phi, float m)
                                             {
                        TLorentzVector p4;
                        p4.SetPtEtaPhiM(pt, eta, phi, m);
                        return p4; }, capturedVar);
                }
                else if (operation == "TLVPtEtaPhiM_corr")
                {
                    // parse the TLorentzVector input
                    auto capturedVar = (this->extractTLVComp)(std::get<2>(step));
                    origRDF = origRDF.Define(std::get<1>(step), [](double pt, float eta, float phi, float m)
                                             {
                        TLorentzVector p4;
                        p4.SetPtEtaPhiM(pt, eta, phi, m);
                        return p4; }, capturedVar);
                }
                else if (operation == "TLVPtEtaPhiE")
                {
                    // parse the TLorentzVector input
                    auto capturedVar = this->extractTLVComp(std::get<2>(step));
                    origRDF = origRDF.Define(std::get<1>(step), [](float pt, float eta, float phi, float e)
                                             {
                        TLorentzVector p4;
                        p4.SetPtEtaPhiE(pt, eta, phi, e);
                        return p4; }, capturedVar);
                }
                else if (operation == "defineDR")
                {
                    // parse the TLorentzVector input
                    auto capturedVar = this->extractTLVComp(std::get<2>(step));
                    origRDF = origRDF.Define(std::get<1>(step), [](float eta1, float eta2, float phi1, float phi2)
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