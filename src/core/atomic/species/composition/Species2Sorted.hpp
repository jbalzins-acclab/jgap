#ifndef JGAP_SPECIES2SORTED_HPP
#define JGAP_SPECIES2SORTED_HPP

#include <algorithm>
#include <array>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <ranges>
#include "core/atomic/species/Species.hpp"
#include "core/atomic/neighbours/NeighbourLists.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {
    struct Species2Sorted {
        static constexpr size_t N = 2;
        std::array<Species, 2> nodes;

        Species2Sorted(const Species2Sorted &) = default;

        Species2Sorted(Species s1, Species s2) : nodes{s1, s2} {
            if (nodes[1] < nodes[0]) std::swap(nodes[1], nodes[0]);
        }

        explicit Species2Sorted(const std::string &encoded) : nodes{Species::Anon(), Species::Anon()} {
            std::stringstream ss(encoded);
            std::string token;
            size_t i = 0;
            while (std::getline(ss, token, ',')) {
                if (i >= 2) JGAP_LOG_AND_THROW("Species2Sorted '{}' has more than the expected 2 node(s)", encoded);
                nodes[i++] = Species(token);
            }
            if (i != 2) JGAP_LOG_AND_THROW("Species2Sorted '{}' has {} node(s), expected 2", encoded, i);
            if (nodes[1] < nodes[0]) std::swap(nodes[1], nodes[0]);
        }

        bool operator==(const Species2Sorted &other) const { return nodes == other.nodes; }
        bool operator<(const Species2Sorted &other) const { return nodes < other.nodes; }

        std::string toString() const {
            return std::string(nodes[0].symbol()) + "," + std::string(nodes[1].symbol());
        }

        static std::set<Species2Sorted> getAll(const NeighbourLists &nl) {
            std::set<Species2Sorted> result;
            for (const auto &[s1, indexes]: nl.atoms_by_species) {
                for (const auto &atom1_idx: indexes) {
                    for (const Species s2: nl.neighbours_per_atom[atom1_idx] | std::views::keys) {
                        result.emplace(s1, s2);
                    }
                }
            }
            return result;
        }
    };
}

#endif
