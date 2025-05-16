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
        {"Run3Summer23NanoAODv12", "Collisions2023_366403_369802_eraBC_GoldenJson"},
        {"Run3Summer23BPixNanoAODv12", "Collisions2023_369803_370790_eraD_GoldenJson"},
    };
    auto it = names.find(era);
    if (it == names.end()) {
        throw std::runtime_error("Unknown era: " + era);
    }
    corr = cset->at(it->second);
}

float PUReweightFunc(float Pileup_nTrueInt) {
    float safe_PU = std::clamp(Pileup_nTrueInt, 0.0f, 99.0f);
    float PUWeight = corr->evaluate({safe_PU, "nominal"});
    return PUWeight;
}