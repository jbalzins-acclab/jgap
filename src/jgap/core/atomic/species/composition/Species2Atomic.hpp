#ifndef JGAP_SPECIES2ATOMIC_HPP
#define JGAP_SPECIES2ATOMIC_HPP

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
    struct Species2Atomic {
        static constexpr size_t N = 2;
        Species root;
        Species node;

        Species2Atomic(const Species2Atomic &) = default;

        Species2Atomic(Species root, Species node) : root(root), node(node) {}

        explicit Species2Atomic(const std::string &encoded) : root(Species::Anon()), node(Species::Anon()) {
            const auto bar = encoded.find('|');
            if (bar == std::string::npos) {
                JGAP_LOG_AND_THROW("Species2Atomic '{}' is missing the '|' root separator", encoded);
            }
            root = Species(encoded.substr(0, bar));
            std::string nodes_part = encoded.substr(bar + 1);
            
            if (nodes_part.empty() || nodes_part.find(',') != std::string::npos) {
                JGAP_LOG_AND_THROW("Species2Atomic '{}' expected exactly 1 node", encoded);
            }
            node = Species(nodes_part);
        }

        bool operator==(const Species2Atomic &other) const {
            return std::tie(root, node) == std::tie(other.root, other.node);
        }
        bool operator<(const Species2Atomic &other) const {
            return std::tie(root, node) < std::tie(other.root, other.node);
        }

        std::string toString() const {
            return std::string(root.symbol()) + "|" + std::string(node.symbol());
        }

        static std::set<Species2Atomic> getAll(const NeighbourLists &nl) {
            std::set<Species2Atomic> result;
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
