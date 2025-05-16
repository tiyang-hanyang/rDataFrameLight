#include <functional>
#include <stdexcept>
#include <string>
#include <map>
#include <memory>
#include <algorithm>

std::shared_ptr<const correction::Correction> corr;

void bind_correction(const std::string& era) { 
    if (!cset) {
        throw std::runtime_error("CorrectionSet not loaded before binding!");
    }
    static const std::map<std::string, std::string> names = {
        {"Run3Summer22EENanoAODv12", "Summer22EE_23Sep2023_RunEFG_V1"},
        {"Run3Summer23NanoAODv12", "Summer23Prompt23_RunC_V1"},
        {"Run2023C", "Summer23Prompt23_RunC_V1"},
        {"Run3Summer23BPixNanoAODv12", "Summer23BPixPrompt23_RunD_V1"},
        {"Run2023D", "Summer23BPixPrompt23_RunD_V1"}
    };
    auto it = names.find(era);
    if (it == names.end()) {
        throw std::runtime_error("Unknown era: " + era);
    }
    corr = cset->at(it->second);
}

ROOT::VecOps::RVec<int> passJetVetoFunc(ROOT::VecOps::RVec<float> eta, ROOT::VecOps::RVec<float> phi) {
    ROOT::VecOps::RVec<int> pass_flags(eta.size());
    for (size_t i = 0; i < eta.size(); ++i) {
        float safe_phi = std::clamp(phi[i], -3.1415f, 3.1415f);
        pass_flags[i] = (corr->evaluate(std::vector<correction::Variable::Type>{"jetvetomap", eta[i], safe_phi}) == 0);
    }
    return pass_flags;
}