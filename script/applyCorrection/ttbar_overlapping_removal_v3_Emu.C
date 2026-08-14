#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <vector>

int TraceOriginV3Emu(const ROOT::VecOps::RVec<Int_t>& pdgId,
                     const ROOT::VecOps::RVec<Short_t>& motherIdx,
                     int targetIdx)
{
    int curIdx = targetIdx;
    const auto nGenPart = pdgId.size();
    const int maxStep = 100;
    for (int step = 0; step < maxStep; ++step)
    {
        if (curIdx < 0 || curIdx >= nGenPart) break;
        const int curPdg = std::abs(pdgId[curIdx]);
        if (curPdg == 6 || curPdg == 25 || curPdg == 23 || curPdg == 24) return curPdg;
        const int mother = motherIdx[curIdx];
        if (mother < 0) break;
        curIdx = mother;
    }
    return 1;
}

bool IsBHadronV3Emu(int pdgId)
{
    const int ap = std::abs(pdgId);
    if (ap < 500) return false;
    const int d3 = (ap / 100) % 10;
    const int d4 = (ap / 1000) % 10;
    return (d3 == 5 || d4 == 5);
}

bool IsLastCopyV3Emu(int statusFlags)
{
    return (statusFlags & (1 << 13)) != 0;
}

ROOT::VecOps::RVec<int> BuildLeptonForOverlapMaskV3Emu(const ROOT::VecOps::RVec<Int_t>& GenPart_pdgId,
                                                       const ROOT::VecOps::RVec<Int_t>& GenPart_status,
                                                       const ROOT::VecOps::RVec<UShort_t>& GenPart_statusFlags,
                                                       const ROOT::VecOps::RVec<Float_t>& GenPart_pt,
                                                       const ROOT::VecOps::RVec<Float_t>& GenPart_eta)
{
    const auto nGenPart = GenPart_pdgId.size();
    ROOT::VecOps::RVec<int> out(nGenPart, 0);
    for (size_t i = 0; i < nGenPart; ++i)
    {
        const int absPdgId = std::abs(GenPart_pdgId[i]);
        const bool isMuon = absPdgId == 13;
        const bool isElectron = absPdgId == 11;
        const bool isStable = GenPart_status[i] == 1;
        const bool isLastCopy = IsLastCopyV3Emu(GenPart_statusFlags[i]);
        const bool passMuonKin = isMuon && (GenPart_pt[i] >= 25.0f) && (std::abs(GenPart_eta[i]) < 2.4f);
        const bool passElectronKin = isElectron && (GenPart_pt[i] >= 30.0f) && (std::abs(GenPart_eta[i]) < 2.5f);
        if ((passMuonKin || passElectronKin) && isStable && isLastCopy) out[i] = 1;
    }
    return out;
}

ROOT::VecOps::RVec<int> GetBHadSourceV3Emu(const ROOT::VecOps::RVec<Int_t>& pdgId,
                                           const ROOT::VecOps::RVec<Short_t>& motherIdx)
{
    const auto nGenPart = pdgId.size();
    std::vector<int> bHadIdx;
    std::unordered_set<int> bHadMotherIdx;
    std::vector<int> terminalBHadIdx;

    for (size_t i = 0; i < nGenPart; ++i)
    {
        if (!IsBHadronV3Emu(pdgId[i])) continue;
        bHadIdx.push_back(i);

        const int mother = motherIdx[i];
        if (mother >= 0 && mother < nGenPart && IsBHadronV3Emu(pdgId[mother])) bHadMotherIdx.insert(mother);
    }

    for (const int idx : bHadIdx)
    {
        if (bHadMotherIdx.find(idx) == bHadMotherIdx.end()) terminalBHadIdx.push_back(idx);
    }

    ROOT::VecOps::RVec<int> origin(nGenPart, 0);
    for (const int idxToSrc : terminalBHadIdx)
    {
        origin[idxToSrc] = TraceOriginV3Emu(pdgId, motherIdx, idxToSrc);
    }
    return origin;
}

float DeltaPhiV3Emu(float phi1, float phi2)
{
    float dPhi = std::abs(phi1 - phi2);
    if (dPhi > 3.1415926f) dPhi = 2.0f * 3.1415926f - dPhi;
    return dPhi;
}

float DeltaRV3Emu(float eta1, float phi1, float eta2, float phi2)
{
    const float dEta = eta1 - eta2;
    const float dPhi = DeltaPhiV3Emu(phi1, phi2);
    return std::sqrt(dEta * dEta + dPhi * dPhi);
}

ROOT::VecOps::RVec<float> GetGenJetMinLeptonDrV3Emu(const ROOT::VecOps::RVec<Float_t>& GenJet_eta,
                                                    const ROOT::VecOps::RVec<Float_t>& GenJet_phi,
                                                    const ROOT::VecOps::RVec<int>& GenPart_isLeptonForOverlap,
                                                    const ROOT::VecOps::RVec<Float_t>& GenPart_eta,
                                                    const ROOT::VecOps::RVec<Float_t>& GenPart_phi)
{
    const auto nGenJet = GenJet_eta.size();
    const auto nGenPart = GenPart_eta.size();
    ROOT::VecOps::RVec<float> out(nGenJet, 999.0f);

    for (size_t jetIdx = 0; jetIdx < nGenJet; ++jetIdx)
    {
        float minDr = 999.0f;
        for (size_t partIdx = 0; partIdx < nGenPart; ++partIdx)
        {
            if (GenPart_isLeptonForOverlap[partIdx] == 0) continue;
            const float dr = DeltaRV3Emu(GenJet_eta[jetIdx], GenJet_phi[jetIdx], GenPart_eta[partIdx], GenPart_phi[partIdx]);
            if (dr < minDr) minDr = dr;
        }
        out[jetIdx] = minDr;
    }
    return out;
}

ROOT::VecOps::RVec<int> GetGenJetBHadV3Emu(const ROOT::VecOps::RVec<UChar_t>& GenJet_hadronFlavour,
                                           const ROOT::VecOps::RVec<Float_t>& GenJet_pt,
                                           const ROOT::VecOps::RVec<Float_t>& GenJet_eta,
                                           const ROOT::VecOps::RVec<Float_t>& GenJet_phi,
                                           const ROOT::VecOps::RVec<int>& bHadron_origin,
                                           const ROOT::VecOps::RVec<Float_t>& GenPart_eta,
                                           const ROOT::VecOps::RVec<Float_t>& GenPart_phi,
                                           const ROOT::VecOps::RVec<float>& GenJet_minLeptonDr,
                                           float leptonDrThreshold = 0.4f)
{
    const auto nGenJet = GenJet_hadronFlavour.size();
    const auto nGenPart = GenPart_eta.size();
    ROOT::VecOps::RVec<int> out(nGenJet, -1);

    for (size_t jetIdx = 0; jetIdx < nGenJet; ++jetIdx)
    {
        if (GenJet_hadronFlavour[jetIdx] != 5) continue;
        if (GenJet_pt[jetIdx] <= 20.0f) continue;
        if (std::abs(GenJet_eta[jetIdx]) >= 2.5f) continue;

        if (GenJet_minLeptonDr[jetIdx] < leptonDrThreshold)
        {
            out[jetIdx] = -13;
            continue;
        }

        int genJetHadID = -2;
        float minDR2 = 0.16f;
        for (size_t partIdx = 0; partIdx < nGenPart; ++partIdx)
        {
            const int partSrc = bHadron_origin[partIdx];
            if (partSrc == 0) continue;

            const float dEta = GenPart_eta[partIdx] - GenJet_eta[jetIdx];
            const float dPhi = DeltaPhiV3Emu(GenPart_phi[partIdx], GenJet_phi[jetIdx]);
            const float dR2 = dEta * dEta + dPhi * dPhi;
            if (dR2 < minDR2)
            {
                minDR2 = dR2;
                genJetHadID = static_cast<int>(partIdx);
            }
        }
        out[jetIdx] = genJetHadID;
    }
    return out;
}

ROOT::VecOps::RVec<int> GetBGenJetsSourceV3Emu(const ROOT::VecOps::RVec<int>& GenJet_hadronIdx,
                                               const ROOT::VecOps::RVec<int>& GenPart_bHadOrigin)
{
    const auto nGenJets = GenJet_hadronIdx.size();
    ROOT::VecOps::RVec<int> out(nGenJets, -1);

    for (size_t genJetIdx = 0; genJetIdx < nGenJets; ++genJetIdx)
    {
        const int bHadIdx = GenJet_hadronIdx[genJetIdx];
        if (bHadIdx == -1 || bHadIdx == -2 || bHadIdx == -13)
        {
            out[genJetIdx] = bHadIdx;
            continue;
        }
        out[genJetIdx] = GenPart_bHadOrigin[bHadIdx];
    }
    return out;
}
