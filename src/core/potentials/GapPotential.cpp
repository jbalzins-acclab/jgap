#include "core/potentials/GapPotential.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    GapPotential::GapPotential(const nlohmann::json &params) {
        CurrentLogger::get()->debug("Parsing jGAP potential params");
        _descriptors = {};
        for (const auto& [label, descriptorParams]: params["descriptors"].items()) {
            _descriptors[label] = ParserRegistry<Descriptor>::get(descriptorParams);
        }
    }

    Predictions GapPotential::predict(const AtomicStructure &structure) {
        Predictions prediction{};
        for (const auto &descriptor: _descriptors | views::values) {
            prediction = prediction + descriptor->predict(structure);
        }
        return prediction;
    }

    nlohmann::json GapPotential::serialize() {
        nlohmann::json descriptors;

        for (const auto &[descriptorLabel, descriptor] : _descriptors) {
            descriptors[descriptorLabel] = descriptor->serialize();
            descriptors[descriptorLabel]["type"] = descriptor->getType();
        }

        return{
            {"descriptors", descriptors}
        };
    }

    CutoffRanges GapPotential::getCutoff() {
        CutoffRanges cutoff{};
        for (const auto& descriptor : _descriptors | views::values) {
            cutoff += descriptor->getCutoff();
        }
        return cutoff;
    }

    void GapPotential::tabulate(TabulationData &table) {
        for (const auto& [label, descriptor]: _descriptors) {
            CurrentLogger::get()->debug("Tabulating {} descriptor", label);
            descriptor->tabulate(table);
        }
    }
}
