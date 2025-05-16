#include <functional>
#include <stdexcept>
#include <string>
#include <map>
#include <memory>
#include <algorithm>

ROOT::VecOps::RVec<float> MuonIDScale(ROOT::VecOps::RVec<float> eta, ROOT::VecOps::RVec<float> pt)
{
    ROOT::VecOps::RVec<float> MuonIDSF(eta.size());
    for (size_t i = 0; i < eta.size(); ++i)
    {
        float safePt = std::clamp(pt[i], 15.1f, 1000.0f);
        float safeEta = std::clamp(eta[i], -2.39f, 2.39f);
        // float safe_eta = std::clamp(eta[i], -2.4f, 2.4f);
        MuonIDSF[i] = cset->at("NUM_TightID_DEN_TrackerMuons")->evaluate(std::vector<correction::Variable::Type>{safeEta, safePt, "nominal"});
    }
    return MuonIDSF;
}

ROOT::VecOps::RVec<float> MuonIsoScale(ROOT::VecOps::RVec<float> eta, ROOT::VecOps::RVec<float> pt)
{
    ROOT::VecOps::RVec<float> MuonIsoSF(eta.size());
    for (size_t i = 0; i < eta.size(); ++i)
    {
        float safePt = std::clamp(pt[i], 15.1f, 1000.0f);
        float safeEta = std::clamp(eta[i], -2.39f, 2.39f);
        // float safe_eta = std::clamp(eta[i], -2.4f, 2.4f);
        MuonIsoSF[i] = cset->at("NUM_TightPFIso_DEN_TightID")->evaluate(std::vector<correction::Variable::Type>{safeEta, safePt, "nominal"});
    }
    return MuonIsoSF;
}


float EventMuonSF(const ROOT::VecOps::RVec<float>& IDSF,
                  const ROOT::VecOps::RVec<float>& IsoSF,
                  const ROOT::VecOps::RVec<int>& good25,
                  const ROOT::VecOps::RVec<int>& good15)
{
    int leading_idx = -1;
    int subleading_idx = -1;

    for (size_t i = 0; i < IDSF.size(); ++i)
    {
        if (good25[i]) {
            leading_idx = i;
            break;
        }
    }

    for (size_t i = 0; i < IDSF.size(); ++i)
    {
        if (good15[i] && i != leading_idx) {
            subleading_idx = i;
            break;
        }
    }

    if (leading_idx < 0 || subleading_idx < 0) return 1.0;
    return IDSF[leading_idx] * IDSF[subleading_idx] * IsoSF[leading_idx] * IsoSF[subleading_idx];
}