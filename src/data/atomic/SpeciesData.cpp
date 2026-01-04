#include "SpeciesData.hpp"

namespace jgap {
    EncodedSpecies SpeciesEncoder::encode(const Species& species) {
        std::lock_guard lock(mtx_);
        if (species_encoder_.count(species)) {
            return species_encoder_.at(species);
        }
        EncodedSpecies encoded = species_counter_++;
        species_encoder_[species] = encoded;
        species_decoder_[encoded] = species;
        return encoded;
    }

    EncodedSpeciesSets SpeciesEncoder::encode(const SpeciesSets &species_sets) {
        EncodedSpeciesSets result{};
        size_t index = 0;
        for (const auto& set: species_sets) {
            assert(index < result.size() && "Species sets too large to be encoded");
            result[index++] = static_cast<std::uint16_t>(set.size());

            for (const auto& species: set) {
                assert(index < result.size() && "Species sets too large to be encoded");
                result[index++] = encode(species);
            }
        }
        // zero = exit
        return result;
    }

    Species SpeciesEncoder::decode(const EncodedSpecies &encoded) {
        std::lock_guard lock(mtx_);
        assert(species_decoder_.contains(encoded) && "Encoded species not found");
        return species_decoder_.at(encoded);
    }

    SpeciesSets SpeciesEncoder::decode(const EncodedSpeciesSets &encoded) {
        SpeciesSets result;
        size_t index = 0;
        while (index < encoded.size() && encoded[index] != 0) { // 0 indicates end of encoded data
            size_t set_size = encoded[index++];
            std::set<Species> current_set;
            for (size_t i = 0; i < set_size; ++i) {
                assert(index < encoded.size()
                    && "Malformed EncodedSpeciesSets: ran out of data while decoding species.");
                current_set.insert(decode(encoded[index++]));
            }
            result.push_back(current_set);
        }
        return result;
    }

    DataNode SpeciesEncoder::toDataNode(const SpeciesSets &species_sets) {
        std::vector<DataNode> result;
        result.reserve(species_sets.size());

        for (const auto& set: species_sets) {
            std::vector<DataNode> species_nodes;
            species_nodes.reserve(set.size());
            for (const auto& species : set) {
                species_nodes.emplace_back(species);
            }
            result.emplace_back(species_nodes);
        }

        return result;
    }

    SpeciesSets SpeciesEncoder::fromDataNode(const DataNode &node) {
        SpeciesSets result;
        if (node.type != DataNode::Type::ARRAY) {
            JGAP_LOG_AND_THROW("DataNode for SpeciesSets must be an array.");
        }

        for (const auto& set_node : node.asArray()) {
            if (set_node.type != DataNode::Type::ARRAY) {
                JGAP_LOG_AND_THROW("Each set in SpeciesSets DataNode must be an array.");
            }
            std::set<Species> current_set;
            for (const auto& species_node : set_node.asArray()) {
                if (species_node.type != DataNode::Type::STRING) {
                    JGAP_LOG_AND_THROW("Each species in a set must be a string.");
                }
                current_set.insert(species_node.asString());
            }
            result.push_back(current_set);
        }
        return result;
    }

    EncodedSpeciesSets SpeciesEncoder::asSet(EncodedSpecies species) {
        EncodedSpeciesSets result{};
        result[0] = 1;
        result[1] = species;
        return result;
    }

    EncodedSpeciesSets SpeciesEncoder::asSet(EncodedSpecies invariant1, EncodedSpecies invariant2) {
        EncodedSpeciesSets result{};
        result[0] = 2;
        result[1] = invariant1;
        result[2] = invariant2;
        return result;
    }

    EncodedSpeciesSets SpeciesEncoder::asSet(EncodedSpecies root, EncodedSpecies node1, EncodedSpecies node2) {
        EncodedSpeciesSets result{};

        result[0] = 1;
        result[1] = root;

        result[2] = 2;
        result[3] = node1;
        result[4] = node2;

        return result;
    }

    double SpeciesEncoder::symmetryFactor(EncodedSpeciesSets &sets) {
        double result = 1.0;
        size_t index = 0;
        while (index < sets.size() && sets[index] != 0) {
            result *= factorial(sets[index]);
            index += sets[index];
        }
        return result;
    }
}
