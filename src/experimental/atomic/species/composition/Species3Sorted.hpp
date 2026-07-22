#ifndef JGAP_SPECIES3SORTED_HPP
#define JGAP_SPECIES3SORTED_HPP

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
    struct Species3Sorted {
        static constexpr size_t N = 3;
        std::array<Species, 3> nodes;

        Species3Sorted(const Species3Sorted &) = default;

        Species3Sorted(Species s1, Species s2, Species s3) : nodes{s1, s2, s3} {
            std::sort(nodes.begin(), nodes.end());
        }

        explicit Species3Sorted(const std::string &encoded) : nodes{Species::Anon(), Species::Anon(), Species::Anon()} {
            std::stringstream ss(encoded);
            std::string token;
            size_t i = 0;
            while (std::getline(ss, token, ',')) {
                if (i >= 3) JGAP_LOG_AND_THROW("Species3Sorted '{}' has more than the expected 3 node(s)", encoded);
                nodes[i++] = Species(token);
            }
            if (i != 3) JGAP_LOG_AND_THROW("Species3Sorted '{}' has {} node(s), expected 3", encoded, i);
            std::sort(nodes.begin(), nodes.end());
        }

        bool operator==(const Species3Sorted &other) const { return nodes == other.nodes; }
        bool operator<(const Species3Sorted &other) const { return nodes < other.nodes; }

        std::string toString() const {
            return std::string(nodes[0].symbol()) + "," + std::string(nodes[1].symbol()) + "," + std::string(nodes[2].symbol());
        }

        static std::set<Species3Sorted> getAll(const NeighbourLists &nl) {
            std::set<Species3Sorted> result;
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
