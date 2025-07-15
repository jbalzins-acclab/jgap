#include "core/potentials/JgapPotential.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    JgapPotential::JgapPotential(const nlohmann::json &params) {
        _descriptors = {};
        for (const auto& [label, descriptorParams]: params["descriptors"].items()) {
            _descriptors[label] = ParserRegistry<Descriptor>::get(descriptorParams);
        }
    }

    PotentialPrediction JgapPotential::predict(const AtomicStructure &structure) {
        PotentialPrediction prediction{};
        for (const auto &descriptor: _descriptors | views::values) {
            prediction = prediction + descriptor->predict(structure);
        }
        return prediction;
    }

    nlohmann::json JgapPotential::serialize() {
        nlohmann::json descriptors;

        for (const auto &[descriptorLabel, descriptor] : _descriptors) {
            descriptors[descriptorLabel] = descriptor->serialize();
            descriptors[descriptorLabel]["type"] = descriptor->getType();
        }

        return{
            {"descriptors", descriptors}
        };
    }

    double JgapPotential::getCutoff() {
        double cutoff = 0.0;
        for (const auto& descriptor : _descriptors | views::values) {
            cutoff = max(cutoff, descriptor->getCutoff());
        }
        return cutoff;
    }

    TabulationData JgapPotential::tabulate(const TabulationParams &params) {
        TabulationData result{};

        for (const auto& descriptor: _descriptors | views::values) {
            result = result + descriptor->tabulate(params);
        }

        return result;
    }
}
