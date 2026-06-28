#ifndef JGAP_SPECIESSET_HPP
#define JGAP_SPECIESSET_HPP

#include "Species.hpp"
#include <array>
#include <algorithm>
#include <string>
#include <sstream>
#include <vector>
#include <initializer_list>

#include "core/atomic/geometry/ClusterSymmetry.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {

    /// \brief Cluster's species stored as [the root if present] + a set/an indexed set of nodes.
    ///
    /// Nodes are sorted upon construction and are immutable afterward,
    /// so given atoms in different relative indexing,
    /// the same cluster would correspond to the same instance of SpeciesSet.
    ///
    /// @note Sorting implies that same-species nodes would have adjacent indices,
    /// which should be exploited when iterating.
    ///
    /// @note Ordering uses internal ids, so order-sensitive data, like coefficients,
    /// should be stored within a node when species sets are used as,
    /// e.g., map keys for these nodes.
    ///
    /// @note The intended workflow is:
    /// define species set of interest -> find all clusters with atom[index].species = set[index]
    /// (where set[0] = root, if root is present), as in \ref NeighbourList.
    /// Finding a cluster first and determining its species set after
    /// requires special care when dealing with indices.
    ///
    /// @tparam NSpecies size of the associated cluster
    /// @tparam ClusterSym permutation symmetry of the cluster
    template<size_t NSpecies, ClusterSymmetry ClusterSym>
    requires(NSpecies > 1)
    class SpeciesSet {
    public:
        static constexpr size_t N = NSpecies;

        // Compile-time determination of layout
        static constexpr bool HasRoot = (ClusterSym != FullSymmetry);
        static constexpr size_t Nodes = NSpecies - HasRoot;

        SpeciesSet(const SpeciesSet& other) = default;

        // Constructor for FullSymmetry (no root, NSpecies nodes)
        template<typename... Args>
        requires(ClusterSym == FullSymmetry && sizeof...(Args) == Nodes)
        explicit SpeciesSet(Args&&... args)
            : nodes{std::forward<Args>(args)...} {
            initNodes();
        }

        // Constructor for NodeSymmetry / IndexSensitive (1 root, NSpecies - 1 nodes)
        template<typename... Args>
        requires(ClusterSym != FullSymmetry && sizeof...(Args) == Nodes)
        SpeciesSet(const Species& r, Args&&... args)
            : root{r}, nodes{std::forward<Args>(args)...} {
            initNodes();
        }

        /// Parses the encoded form produced by \ref SpeciesSet::toString:
        /// FullSymmetry as "n0,n1,...", otherwise "root|n0,n1,...".
        explicit SpeciesSet(const std::string& encoded) {
            std::string nodes_part = encoded;

            if constexpr (HasRoot) {
                const auto bar = encoded.find('|');
                if (bar == std::string::npos) {
                    JGAP_LOG_AND_THROW("SpeciesSet '{}' is missing the '|' root separator", encoded);
                }
                root[0] = Species(encoded.substr(0, bar));
                nodes_part = encoded.substr(bar + 1);
            }

            std::stringstream ss(nodes_part);
            std::string token;
            size_t i = 0;
            while (std::getline(ss, token, ',')) {
                if (i >= Nodes) {
                    JGAP_LOG_AND_THROW("SpeciesSet '{}' has more than the expected {} node(s)", encoded, Nodes);
                }
                nodes[i++] = Species(token);
            }
            if (i != Nodes) {
                JGAP_LOG_AND_THROW("SpeciesSet '{}' has {} node(s), expected {}", encoded, i, Nodes);
            }

            initNodes();
        }

        bool operator==(const SpeciesSet& other) const {
            return std::tie(root, nodes) == std::tie(other.root, other.nodes);
        }

        bool operator<(const SpeciesSet& other) const {
            return std::tie(root, nodes) < std::tie(other.root, other.nodes);
        }

        Species getRoot() const requires(HasRoot) {
            return root[0];
        }

        const std::array<Species, Nodes>& getNodes() const {
            return nodes;
        }

        std::string toString() const {
            std::stringstream ss;
            if constexpr (HasRoot) {
                ss << root[0].symbol() << "|";
            }
            for (size_t i = 0; i < Nodes; i++) {
                ss << (i > 0 ? "," : "") << nodes[i].symbol();
            }
            return ss.str();
        }

    private:
        void initNodes() {

            if constexpr (ClusterSym == IndexSensitive) {
                return;
            }

            if constexpr (Nodes == 2) {
                if (nodes[1] < nodes[0]) std::swap(nodes[1], nodes[0]);
            } else if constexpr (Nodes > 2) {
                std::sort(nodes.begin(), nodes.end());
            }
        }

        [[no_unique_address]] std::array<Species, HasRoot> root{};
        std::array<Species, Nodes> nodes;
    };
}

#endif