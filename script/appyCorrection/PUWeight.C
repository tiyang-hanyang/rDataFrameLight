#include <functional>
#include <stdexcept>
#include <string>
#include <map>
#include <memory>
#include <algorithm>

std::shared_ptr<const correction::CorrectionSet> PU_cset;
std::shared_ptr<const correction::Correction> PU_corr;

void PU_weight_init(const std::string& era) { 
    // hard coding the path of the file here (as just macro)
    static const std::map<std::string, std::string> csetFile = {
        {"RunIII2024Summer24NanoAODv15", "/home/tiyang/public/rDataFrameLight_git/correction/LUM/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/puWeights_BCDEFGHI.json"},
        {"Run3Summer23BPixNanoAODv12", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/LUM/2023_Summer23BPix/puWeights.json"},
        {"Run3Summer23NanoAODv12", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/LUM/2023_Summer23/puWeights.json"},
        {"Run3Summer22NanoAODv12", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/LUM/2022_Summer22/puWeights.json"},
        {"Run3Summer22EENanoAODv12", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/LUM/2022_Summer22EE/puWeights.json"},
    };

    PU_cset = correction::CorrectionSet::from_file(csetFile.at(era));

    if (!PU_cset) {
        throw std::runtime_error("CorrectionSet not loaded before binding!");
    }
    static const std::map<std::string, std::string> names = {
        {"Run3Summer23NanoAODv12", "Collisions2023_366403_369802_eraBC_GoldenJson"},
        {"Run3Summer23BPixNanoAODv12", "Collisions2023_369803_370790_eraD_GoldenJson"},
        {"Run3Summer22NanoAODv12", "Collisions2022_355100_357900_eraBCD_GoldenJson"},
        {"Run3Summer22EENanoAODv12", "Collisions2022_359022_362760_eraEFG_GoldenJson"},
        {"RunIII2024Summer24NanoAODv15", "Collisions24_BCDEFGHI_goldenJSON"}
        // or if only C
        // {"RunIII2024Summer24NanoAODv15", "Collisions24_C_goldenJSON"}
    };
    auto it = names.find(era);
    if (it == names.end()) {
        throw std::runtime_error("Unknown era: " + era);
    }
    PU_corr = PU_cset->at(it->second);
}

float PUReweightFunc(float Pileup_nTrueInt) {
    float safe_PU = std::clamp(Pileup_nTrueInt, 0.0f, 99.0f);
    float PUWeight = PU_corr->evaluate({safe_PU, "nominal"});
    return PUWeight;
}