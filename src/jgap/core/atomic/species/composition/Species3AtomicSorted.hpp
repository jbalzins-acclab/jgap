#ifndef JGAP_SPECIES3ATOMICSORTED_HPP
#define JGAP_SPECIES3ATOMICSORTED_HPP

#include <algorithm>
#include <array>
#include <tuple>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <ranges>
#include "jgap/core/atomic/species/Species.hpp"
#include "jgap/core/atomic/neighbours/NeighbourLists.hpp"
#include "../../../io/log/CurrentLogger.hpp"

namespace jgap {
    struct Species3AtomicSorted {
        static constexpr size_t N = 3;
        Species root;
        std::array<Species, 2> nodes;

        Species3AtomicSorted(const Species3AtomicSorted &) = default;

        Species3AtomicSorted(Species root, Species s1, Species s2) : root(root), nodes{s1, s2} {
            if (nodes[1] < nodes[0]) std::swap(nodes[1], nodes[0]);
        }

        explicit Species3AtomicSorted(const std::string &encoded) : root(Species::Anon()), nodes{Species::Anon(), Species::Anon()} {
            const auto bar = encoded.find('|');
            if (bar == std::string::npos) {
                JGAP_LOG_AND_THROW("Species3AtomicSorted '{}' is missing the '|' root separator", encoded);
            }
            root = Species(encoded.substr(0, bar));
            std::string nodes_part = encoded.substr(bar + 1);

            std::stringstream ss(nodes_part);
            std::string token;
            size_t i = 0;
            while (std::getline(ss, token, ',')) {
                if (i >= 2) JGAP_LOG_AND_THROW("Species3AtomicSorted '{}' has more than the expected 2 node(s)", encoded);
                nodes[i++] = Species(token);
            }
            if (i != 2) JGAP_LOG_AND_THROW("Species3AtomicSorted '{}' has {} node(s), expected 2", encoded, i);
            if (nodes[1] < nodes[0]) std::swap(nodes[1], nodes[0]);
        }

        bool operator==(const Species3AtomicSorted &other) const {
            return std::tie(root, nodes) == std::tie(other.root, other.nodes);
        }
        bool operator<(const Species3AtomicSorted &other) const {
            return std::tie(root, nodes) < std::tie(other.root, other.nodes);
        }

        std::string toString() const {
            return std::string(root.symbol()) + "|" + std::string(nodes[0].symbol()) + "," + std::string(nodes[1].symbol());
        }

        static std::set<Species3AtomicSorted> getAll(const NeighbourLists &nl) {
            std::set<Species3AtomicSorted> result;
            for (const auto &[s1, indexes]: nl.atoms_by_species) {
                for (const auto &atom1_idx: indexes) {
                    for (const Species s2: nl.neighbours_per_atom[atom1_idx] | std::views::keys) {
                        for (const Species s3: nl.neighbours_per_atom[atom1_idx] | std::views::keys) {
                            result.emplace(s1, s2, s3);
                        }
                    }
                }
            }
            return result;
        }
    };
}

#endif
