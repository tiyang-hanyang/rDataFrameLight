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
        {"Run2022C", "Summer22_23Sep2023_RunCD_V1"},
        {"Run2022D", "Summer22_23Sep2023_RunCD_V1"},
        {"Run3Summer22NanoAODv12", "Summer22_23Sep2023_RunCD_V1"},
        {"Run2022E", "Summer22EE_23Sep2023_RunEFG_V1"},
        {"Run2022F", "Summer22EE_23Sep2023_RunEFG_V1"},
        {"Run2022G", "Summer22EE_23Sep2023_RunEFG_V1"},
        {"Run3Summer22EENanoAODv12", "Summer22EE_23Sep2023_RunEFG_V1"},
        {"Run3Summer23NanoAODv12", "Summer23Prompt23_RunC_V1"},
        {"Run2023C", "Summer23Prompt23_RunC_V1"},
        {"Run3Summer23BPixNanoAODv12", "Summer23BPixPrompt23_RunD_V1"},
        {"Run2023D", "Summer23BPixPrompt23_RunD_V1"},
        {"Run2024C", "Summer24Prompt24_RunBCDEFGHI_V1"},
        {"Run2024D", "Summer24Prompt24_RunBCDEFGHI_V1"},
        {"Run2024E", "Summer24Prompt24_RunBCDEFGHI_V1"},
        {"Run2024F", "Summer24Prompt24_RunBCDEFGHI_V1"},
        {"Run2024G", "Summer24Prompt24_RunBCDEFGHI_V1"},
        {"Run2024H", "Summer24Prompt24_RunBCDEFGHI_V1"},
        {"Run2024I", "Summer24Prompt24_RunBCDEFGHI_V1"},
        {"RunIII2024Summer24NanoAODv15", "Summer24Prompt24_RunBCDEFGHI_V1"}
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
        // exclude non coverage etas
        if (abs(eta[i])>2.5) {
            pass_flags[i]=0;
            continue;
        }
        float safe_phi = std::clamp(phi[i], -3.1415f, 3.1415f);
        pass_flags[i] = (corr->evaluate(std::vector<correction::Variable::Type>{"jetvetomap", eta[i], safe_phi}) == 0);
    }
    return pass_flags;
}