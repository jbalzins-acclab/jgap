#include "core/potentials/JtabGapPotential.hpp"

namespace jgap {
    JtabGapPotential::JtabGapPotential(const nlohmann::json &params) {
    }

    PotentialPrediction JtabGapPotential::predict(const AtomicStructure &structure) {
        throw 0;
    }

    nlohmann::json JtabGapPotential::serialize() {
        throw 0;
    }

    double JtabGapPotential::getCutoff() {
        return 0;
    }

    TabulationData JtabGapPotential::tabulate(const TabulationParams &params) {
        return {};
    }
}
