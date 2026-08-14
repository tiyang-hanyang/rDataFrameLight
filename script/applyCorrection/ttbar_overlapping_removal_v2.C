#include <functional>
#include <stdexcept>
#include <string>
#include <map>
#include <memory>
#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <vector>

//////////////////////////////////////////////////
// assist functions
// for getting the source of all b hadrons
// iterative tracing the decay chain
// record as abs values
// I need to distinguish the case of not b-hadron, or b-hadron not from WZHt
// get into this function cannot return 0
// 0 is reserved for "not a terminal b-hadron"
// use 1 for a matched terminal b-hadron whose origin is neither 6/23/24/25
int TraceOrigin(const ROOT::VecOps::RVec<Int_t>& pdgId,
                const ROOT::VecOps::RVec<Short_t>& motherIdx,
                int targetIdx)
{
    int curIdx = targetIdx;
    auto nGenPart = pdgId.size();
    int maxStep = 100;
    for (int step = 0; step < maxStep; ++step) 
    {
        if (curIdx < 0 || curIdx >= nGenPart) break;
        int curPdg = std::abs(pdgId[curIdx]);
        if (curPdg == 6 || curPdg == 25 || curPdg == 23 || curPdg == 24) return curPdg;
        int mother = motherIdx[curIdx];
        if (mother < 0) break;
        curIdx = mother;
    }
    return 1;
}

// b meson 5xx, b baryon 5xxx
bool IsBHadron(int pdgId) 
{
    int ap = std::abs(pdgId);
    if (ap < 500) return false;
    int d3 = (ap / 100) % 10;
    int d4 = (ap / 1000) % 10;
    return (d3 == 5 || d4 == 5);
}

//////////////////////////////////////////////////
// utilize
// first, need to trace all the terminal b-hadron's source
// 0 means not a terminal b-hadron
// 1 means b-hadron not from 6 23 24 25
ROOT::VecOps::RVec<int> GetBHadSource(const ROOT::VecOps::RVec<Int_t>& pdgId, const ROOT::VecOps::RVec<Short_t>& motherIdx)
{
    auto nGenPart = pdgId.size();
    std::vector<int> bHadIdx;
    std::unordered_set<int> bHadMotherIdx;
    std::vector<int> terminalBHadIdx;

    for (size_t i = 0; i < nGenPart; ++i)
    {
        if (!IsBHadron(pdgId[i])) continue;
        bHadIdx.push_back(i);

        int mother = motherIdx[i];
        if (mother >= 0 && mother < nGenPart && IsBHadron(pdgId[mother]))
        {
            bHadMotherIdx.insert(mother);
        }
    }

    for (size_t i = 0; i < bHadIdx.size(); ++i)
    {
        int tempIdx = bHadIdx[i];
        if (bHadMotherIdx.find(tempIdx) == bHadMotherIdx.end()) terminalBHadIdx.push_back(tempIdx);
    }

    ROOT::VecOps::RVec<int> origin(nGenPart, 0);
    auto nTerminal = terminalBHadIdx.size();
    for (size_t i = 0; i < nTerminal; ++i)
    {
        auto idxToSrc = terminalBHadIdx[i];
        origin[idxToSrc] = TraceOrigin(pdgId, motherIdx, idxToSrc);
    }
    return origin;
}

//////////////////////////////////////////////////
// utilize
// then, let us get the deltaR matching to the GenParts
// -1 means not in the fiducial b-jet pool
//    i.e. not (hadronFlavour==5 and pt>20 and |eta|<2.5)
// -2 means in the fiducial b-jet pool but not matched to any terminal b-hadron
ROOT::VecOps::RVec<int> GetGenJetBHad(const ROOT::VecOps::RVec<UChar_t>& GenJet_hadronFlavour, const ROOT::VecOps::RVec<Float_t>& GenJet_pt, const ROOT::VecOps::RVec<Float_t>& GenJet_eta, const ROOT::VecOps::RVec<Float_t>& GenJet_phi, const ROOT::VecOps::RVec<int>& bHadron_origin, const ROOT::VecOps::RVec<Float_t>& GenPart_eta, const ROOT::VecOps::RVec<Float_t>&GenPart_phi)
{
    auto genJet_size = GenJet_hadronFlavour.size();
    auto nGenPart = GenPart_eta.size();
    // default: not in the fiducial b-jet pool
    ROOT::VecOps::RVec<int> GenJet_hadronIdx(genJet_size, -1);
    for (size_t i = 0; i < genJet_size; ++i) 
    {
        // only match for fiducial bjets
        if ( GenJet_hadronFlavour[i] != 5 ) continue;
        if ( GenJet_pt[i] <= 20.0f ) continue;
        if ( std::abs(GenJet_eta[i]) >= 2.5f ) continue;

        float genJetEta = GenJet_eta[i];
        float genJetPhi = GenJet_phi[i];

        // try to match the genpart with the smallest detlaR, threshold 0.4
        int genJetHadID = -2;
        float minDR2 = 0.16;
        for (size_t partIdx =0; partIdx < nGenPart; ++partIdx)
        {
            int partSrc = bHadron_origin[partIdx];
            // only filter out non-terminal/non-b-hadron (0)
            // keep traced sources 6/23/24/25 and also the generic non-top code (1)
            if (partSrc == 0) continue; 
            float hadEta = GenPart_eta[partIdx];
            float hadPhi = GenPart_phi[partIdx];

            float dEta2 = (hadEta - genJetEta) * (hadEta - genJetEta);
            float dPhi = std::abs(hadPhi - genJetPhi);
            if (dPhi > 3.1415926) dPhi = 2*3.1415926 - dPhi;
            
            float dR2 = dEta2 + dPhi*dPhi;
            if (dR2 < minDR2) 
            {
                minDR2 = dR2;
                genJetHadID = partIdx;
            }
        }
        GenJet_hadronIdx[i] = genJetHadID;
    }
    return GenJet_hadronIdx;
}

ROOT::VecOps::RVec<int> GetGenJetBHad_withoutGM(const ROOT::VecOps::RVec<Float_t>& GenJet_eta, const ROOT::VecOps::RVec<Float_t>& GenJet_phi, const ROOT::VecOps::RVec<int>& bHadron_origin, const ROOT::VecOps::RVec<Float_t>& GenPart_eta, const ROOT::VecOps::RVec<Float_t>&GenPart_phi)
{
    auto genJet_size = GenJet_eta.size();
    auto nGenPart = GenPart_eta.size();
    // the default value is in case if the genjet did not successfully match to any bhadron (dR is too large)
    ROOT::VecOps::RVec<int> GenJet_hadronIdx(genJet_size, -1);
    for (size_t i = 0; i < genJet_size; ++i) 
    {
        // only match for bjets
        // in case if not having this remains
        // if ( GenJet_hadronFlavour[i] != 5 ) continue;

        float genJetEta = GenJet_eta[i];
        float genJetPhi = GenJet_phi[i];

        // try to match the genpart with the smallest detlaR, threshold 0.4
        int genJetHadID = -1;
        float minDR2 = 0.16;
        for (size_t partIdx =0; partIdx < nGenPart; ++partIdx)
        {
            int partSrc = bHadron_origin[partIdx];
            // only filter out non-terminal/non-b-hadron (0)
            // keep traced sources 6/23/24/25 and also the generic non-top code (1)
            if (partSrc == 0) continue; 
            float hadEta = GenPart_eta[partIdx];
            float hadPhi = GenPart_phi[partIdx];

            float dEta2 = (hadEta - genJetEta) * (hadEta - genJetEta);
            float dPhi = std::abs(hadPhi - genJetPhi);
            if (dPhi > 3.1415926) dPhi = 2*3.1415926 - dPhi;
            
            float dR2 = dEta2 + dPhi*dPhi;
            if (dR2 < minDR2) 
            {
                minDR2 = dR2;
                genJetHadID = partIdx;
            }
        }
        GenJet_hadronIdx[i] = genJetHadID;
    }
    return GenJet_hadronIdx;
}


//////////////////////////////////////////////////
// utilize
// finally, try to match the jets to genjets, for identifying the source
// here not needed as the removal should be at gen level
ROOT::VecOps::RVec<int> GetBJetsSource(const ROOT::VecOps::RVec<Short_t>& Jet_genJetIdx, const ROOT::VecOps::RVec<int>& GenJet_hadronIdx, const ROOT::VecOps::RVec<int>& GenPart_bHadOrigin)
{
    auto nJets = Jet_genJetIdx.size();
    auto nGenJets = GenJet_hadronIdx.size();
    ROOT::VecOps::RVec<int> Jets_bsource(nJets, -1);
    // if not in the fiducial b-jet pool, then -1
    // if in the fiducial b-jet pool but not matched to a terminal b-hadron, then -2
    for (size_t jetIdx = 0; jetIdx < nJets; ++jetIdx) 
    {
        int genJetIdx = Jet_genJetIdx[jetIdx];
        if ( genJetIdx<0 || genJetIdx>=nGenJets ) continue;
        int bHadIdx = GenJet_hadronIdx[genJetIdx];
        if (bHadIdx == -1) continue;
        if (bHadIdx == -2)
        {
            Jets_bsource[jetIdx] = -2;
            continue;
        }
        int bHadSrc = GenPart_bHadOrigin[bHadIdx];
        Jets_bsource[jetIdx] = bHadSrc;
    }
    return Jets_bsource;
}

//////////////////////////////////////////////////
// utilize
// finally, try to match the jets to genjets, for identifying the source
// -1 means not in the fiducial b-jet pool
// -2 means fiducial b-jet but not matched to any terminal b-hadron
// 1, 6, 23, 24, 25 are traced sources
ROOT::VecOps::RVec<int> GetBGenJetsSource(const ROOT::VecOps::RVec<int>& GenJet_hadronIdx, const ROOT::VecOps::RVec<int>& GenPart_bHadOrigin)
{
    auto nGenJets = GenJet_hadronIdx.size();
    ROOT::VecOps::RVec<int> GenJets_bsource(nGenJets, -1);

    for (size_t genJetIdx = 0; genJetIdx < nGenJets; genJetIdx++) 
    {
        int bHadIdx = GenJet_hadronIdx[genJetIdx];
        if (bHadIdx == -1) continue;
        if (bHadIdx == -2)
        {
            GenJets_bsource[genJetIdx] = -2;
            continue;
        }
        int bHadSrc = GenPart_bHadOrigin[bHadIdx];
        GenJets_bsource[genJetIdx] = bHadSrc;
    }
    return GenJets_bsource;
}

//////////////////////////////////////////////////
// The problem is the extraction does not really work, as in the ttbar NANOAOD sample, not all information about the decay chain are shown: I would also need to keep one additional method: deviation by the parton level information directly
// In this way, I would just count the multiplicity of the b-quark as the first copy, and look if these are from the top, since in principle this step will always appear inside the decay chain of NANOAOD

bool IsFirstCopy(int statusFlags) 
{
    return statusFlags & (1 << 12);
}

int countAddB(const ROOT::VecOps::RVec<Int_t>& GenPart_pdgId, const ROOT::VecOps::RVec<Short_t>& GenPart_genPartIdxMother, const ROOT::VecOps::RVec<UShort_t>& GenPart_statusFlags)
{
    size_t nPart = GenPart_pdgId.size();
    int nNonTopB(0);
    for (size_t i = 0; i < nPart; i++)
    {
        // check if b-quark
        if (GenPart_pdgId[i] != 5 && GenPart_pdgId[i] != -5) continue;
        // require to be the first particle
        if (!IsFirstCopy(GenPart_statusFlags[i])) continue;

        // to check if the mother particle exist and if that is top quark
        auto motherId = GenPart_genPartIdxMother[i];
        if (motherId == -1) continue;
        auto mother_pdgID = GenPart_pdgId[motherId];
        if (mother_pdgID!=6 && mother_pdgID!=-6) nNonTopB++;
    }
    return nNonTopB;
}
