#ifndef JGAP_SPECIESDATA_HPP
#define JGAP_SPECIESDATA_HPP

#include <algorithm>

#include "Vector3.hpp"

namespace jgap {

    using Species = string;

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

        string toString() const {
            return _pair.first + "," + _pair.second;
        }

        bool contains(const Species& species) const {
            return _pair.first == species || _pair.second == species;
        }

    private:
        pair<Species, Species> _pair;
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

        string toString() const {
            return root + "," + nodes.toString();
        }
    };
}

#endif