#include <algorithm>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

std::shared_ptr<const correction::CorrectionSet> TriggerScale_cset;
std::shared_ptr<const correction::Correction> TriggerScale_corr;
bool TriggerScale_use2D = false;

void TriggerScale_init(const std::string &filePath, const std::string &correctionName, bool use2D = false)
{
    TriggerScale_cset = correction::CorrectionSet::from_file(filePath);
    TriggerScale_corr = TriggerScale_cset->at(correctionName);
    TriggerScale_use2D = use2D;
}

float LeadingMuonTriggerScale(float pt, const std::string &variation = "nominal")
{
    if (!TriggerScale_corr)
    {
        throw std::runtime_error("TriggerScale correction is not initialized.");
    }
    const float safePt = std::max(pt, 0.0f);
    return TriggerScale_corr->evaluate(std::vector<correction::Variable::Type>{safePt, variation});
}

float DiMuonTriggerScale(float leadingPt, float subleadingPt, const std::string &variation = "nominal")
{
    if (!TriggerScale_corr)
    {
        throw std::runtime_error("TriggerScale correction is not initialized.");
    }
    const float safeLeadingPt = std::max(leadingPt, 0.0f);
    const float safeSubleadingPt = std::max(subleadingPt, 0.0f);
    return TriggerScale_corr->evaluate(std::vector<correction::Variable::Type>{safeLeadingPt, safeSubleadingPt, variation});
}

float LeadingMuonTriggerScaleFromRVec(const ROOT::VecOps::RVec<float> &pt, int leadingMuonIdx, const std::string &variation = "nominal")
{
    if (leadingMuonIdx < 0 || static_cast<size_t>(leadingMuonIdx) >= pt.size())
    {
        return 1.0f;
    }
    return LeadingMuonTriggerScale(pt[leadingMuonIdx], variation);
}

float DiMuonTriggerScaleFromRVec(
    const ROOT::VecOps::RVec<float> &pt,
    int leadingMuonIdx,
    int subleadingMuonIdx,
    const std::string &variation = "nominal")
{
    if (leadingMuonIdx < 0 || subleadingMuonIdx < 0)
    {
        return 1.0f;
    }
    if (static_cast<size_t>(leadingMuonIdx) >= pt.size() || static_cast<size_t>(subleadingMuonIdx) >= pt.size())
    {
        return 1.0f;
    }
    return DiMuonTriggerScale(pt[leadingMuonIdx], pt[subleadingMuonIdx], variation);
}

float EventTriggerScaleFromRVec(
    const ROOT::VecOps::RVec<float> &pt,
    int leadingMuonIdx,
    int subleadingMuonIdx,
    const std::string &variation = "nominal")
{
    if (TriggerScale_use2D)
    {
        return DiMuonTriggerScaleFromRVec(pt, leadingMuonIdx, subleadingMuonIdx, variation);
    }
    return LeadingMuonTriggerScaleFromRVec(pt, leadingMuonIdx, variation);
}
