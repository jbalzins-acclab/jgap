#include "core/potentials/GapPotential.hpp"
#include "io/parse/ParserRegistry.hpp"
#include "utils/Utils.hpp"

namespace jgap {
    GapPotential::GapPotential(const DataNode &params) {
        JGAP_LOG_DEBUG("Parsing jGAP potential params");
        _descriptors = {};
        const auto& descsNode = REQUIRE(params, "descriptors");
        const auto& m = std::get<std::map<std::string, DataNode>>(descsNode.value);
        for (const auto& [label, descriptorParams] : m) {
            _descriptors[label] = REGISTRY_GET(Descriptor, descriptorParams);
        }
    }

    Predictions GapPotential::predict(const AtomicStructure &structure) {
        Predictions prediction{};
        for (const auto &descriptor: _descriptors | std::views::values) {
            prediction = prediction + descriptor->predict(structure);
        }
        return prediction;
    }

    DataNode GapPotential::serialize() {
        DataNode descriptors = DataNode::object();
        auto& dm = std::get<std::map<std::string, DataNode>>(descriptors.value);
        for (const auto &[descriptorLabel, descriptor] : _descriptors) {
            descriptors[descriptorLabel] = descriptor->serialize();
            descriptors[descriptorLabel]["type"] = descriptor->getType();
        }

        return DataNode{
            {"descriptors", descriptors}
        };
    }

    CutoffRanges GapPotential::getCutoff() {
        CutoffRanges cutoff{};
        for (const auto& descriptor : _descriptors | std::views::values) {
            cutoff += descriptor->getCutoff();
        }
        return cutoff;
    }

    void GapPotential::tabulate(TabulationData &table) {
        for (const auto& [label, descriptor]: _descriptors) {
            JGAP_LOG_DEBUG("Tabulating {} descriptor", label);
            descriptor->tabulate(table);
        }
    }
}
