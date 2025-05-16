#include <functional>
#include <stdexcept>
#include <string>
#include <map>
#include <memory>
#include <algorithm>

std::shared_ptr<const correction::Correction> corr;

void bind_correction() {
    if (!cset) {
        throw std::runtime_error("CorrectionSet not loaded before binding!");
    }
    corr = cset->at("DY_pTll_reweighting");
}