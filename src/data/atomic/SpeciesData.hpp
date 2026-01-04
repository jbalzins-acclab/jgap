#ifndef JGAP_SPECIESDATA_HPP
#define JGAP_SPECIESDATA_HPP

#include <algorithm>
#include <set>
#include <atomic>
#include <mutex>

#include "../Vector3.hpp"
#include "data/DataNode.hpp"
#include "utils/Utils.hpp"

#define ENCODED_SPECIES_SET_DS_LENGTH 32 // optimized for x86 Cache lane (?)

namespace jgap {

    using Species = std::string;
    using EncodedSpecies = std::uint16_t;

    using SpeciesSets = std::vector<std::set<Species>>;
    using EncodedSpeciesSets = std::array<EncodedSpecies, ENCODED_SPECIES_SET_DS_LENGTH>;

    class SpeciesEncoder {
    public:
        static EncodedSpecies encode(const Species& species);
        static EncodedSpeciesSets encode(const SpeciesSets& species_sets);

        static Species decode(const EncodedSpecies& encoded);
        static SpeciesSets decode(const EncodedSpeciesSets& encoded);

        static DataNode toDataNode(const SpeciesSets& species_sets);
        static SpeciesSets fromDataNode(const DataNode& node);

        static EncodedSpeciesSets asSet(EncodedSpecies species);
        static EncodedSpeciesSets asSet(EncodedSpecies invariant1, EncodedSpecies invariant2);
        static EncodedSpeciesSets asSet(EncodedSpecies root, EncodedSpecies node1, EncodedSpecies node2);

        static double symmetryFactor(EncodedSpeciesSets& sets);

    private:
        inline static std::atomic_uint16_t species_counter_ = 0;
        inline static std::map<Species, EncodedSpecies> species_encoder_{};
        inline static std::map<EncodedSpecies, Species> species_decoder_{};
        inline static std::mutex mtx_{};
    };

    struct ContributorReceiverSpecies {
        Species contributor;
        Species receiver;

        bool operator==(const ContributorReceiverSpecies &other) const {
            return contributor == other.contributor && receiver == other.receiver;
        }

        bool operator<(const ContributorReceiverSpecies &other) const {
            return tie(contributor, receiver) < tie(other.contributor, other.receiver);
        }
    };

    class SpeciesPair {
    public:
        SpeciesPair(const Species& first, const Species& second)
            : _pair(minmax(first, second)) {}

        const Species& first() const { return _pair.first; }
        const Species& second() const { return _pair.second; }

        bool operator==(const SpeciesPair& other) const {
            return _pair == other._pair;
        }

        bool operator<(const SpeciesPair& other) const {
            return _pair < other._pair;
        }

        std::string toString() const {
            return _pair.first + "," + _pair.second;
        }

        bool contains(const Species& species) const {
            return _pair.first == species || _pair.second == species;
        }

    private:
        std::pair<Species, Species> _pair;
    };

    struct SpeciesTriplet {

        Species root;
        SpeciesPair nodes;

        bool operator==(const SpeciesTriplet& other) const {
            return root == other.root && nodes == other.nodes;
        }

        bool operator<(const SpeciesTriplet& other) const {
            if (root == other.root) {
                return nodes < other.nodes;
            }
            return root < other.root;
        }

        std::string toString() const {
            return root + "," + nodes.toString();
        }
    };
}

#endif