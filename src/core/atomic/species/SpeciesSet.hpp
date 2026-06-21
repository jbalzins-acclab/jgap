#ifndef JGAP_SPECIESSET_HPP
#define JGAP_SPECIESSET_HPP

#include "Species.hpp"
#include <array>
#include <algorithm>
#include <string>
#include <sstream>
#include <vector>
#include <initializer_list>
#include "io/log/CurrentLogger.hpp"

namespace jgap {

    /**
     * Indicates whether an atomic cluster({@link Cluster} @ index=0) has a special "root" atom
     * that is not compliant with PERMUTATION symmetry.
     *
     * If yes, this also indicates that the {@link NBodyDescriptor} derived from this cluster
     * is exclusively "owned" by this "root" atom.
     */
    enum class ClusterSymmetry {
        HasCentralAtom,
        FullSymmetry
    };

    using ClusterSymmetry::HasCentralAtom;
    using ClusterSymmetry::FullSymmetry;

    /**
     * Cluster's species stored as [root if present] + set of nodes.
     * Nodes are sorted upon construction and are immutable afterward,
     * so given atoms in different relative indexing,
     * the same cluster would correspond to the same instance of SpeciesSet.
     *
     * @note Ordering uses internal ids, so order-sensitive data, like coefficients,
     * should be stored within a node when species sets are used as,
     * e.g., map keys for these nodes.
     *
     * @note Intended workflow is:
     * define species set of interest -> find all clusters with atom[index].species = set[index]
     * (where set[0] = root, if root is present), as in {@link NeighbourList}.
     * Finding a cluster first and determining its species set after
     * requires special care when dealing with indices.
     *
     * @tparam NSpecies size of the associated cluster
     * @tparam ClusterSym permutation symmetry of the cluster
     */
    template<size_t NSpecies, ClusterSymmetry ClusterSym>
    class SpeciesSet;

    template<size_t NSpecies>
    requires(NSpecies > 1)
    class SpeciesSet<NSpecies, FullSymmetry> {
    public:
        static constexpr size_t N = NSpecies;

        SpeciesSet(const SpeciesSet& other) = default;

        template<typename... Args>
        requires(sizeof...(Args) == NSpecies)
        explicit SpeciesSet(Args&&... args) : nodes{std::forward<Args>(args)...} {
            if constexpr (NSpecies == 2) {
                if (nodes[1] < nodes[0]) {
                    std::swap(nodes[1], nodes[0]);
                }
                return;
            }
            std::sort(nodes.begin(), nodes.end());
        }

        bool operator==(const SpeciesSet& other) const {
            return nodes == other.nodes;
        }

        bool operator<(const SpeciesSet &other) const {
            return nodes < other.nodes;
        }

        const std::array<Species, NSpecies>& getNodes() const {
            return nodes;
        }

        std::string toString() const {
            std::stringstream ss;
            for (size_t i = 0; i < NSpecies; i++) {
                ss << (i > 0 ? "," : "") << nodes[i].symbol();
            }
            return ss.str();
        }

    private:
        std::array<Species, NSpecies> nodes;
    };

    template<size_t NSpecies>
    requires(NSpecies > 1)
    class SpeciesSet<NSpecies, HasCentralAtom> {
    public:
        static constexpr size_t N = NSpecies;

        SpeciesSet(const SpeciesSet& other) = default;

        template<typename... Args>
        requires(sizeof...(Args) == NSpecies - 1)
        SpeciesSet(const Species& root, Args&&... args)
            : root(root), nodes{std::forward<Args>(args)...} {

            if constexpr (NSpecies == 2) {
                return;
            }
            if constexpr (NSpecies == 3) {
                if (nodes[1] < nodes[0]) {
                    std::swap(nodes[1], nodes[0]);
                }
                return;
            }
            std::sort(nodes.begin(), nodes.end());
        }

        bool operator==(const SpeciesSet& other) const {
            return std::tie(root, nodes) == std::tie(other.root, other.nodes);
        }

        bool operator<(const SpeciesSet &other) const {
            return std::tie(root, nodes) < std::tie(other.root, other.nodes);
        }

        const Species& getRoot() const {
            return root;
        }

        const std::array<Species, NSpecies - 1>& getNodes() const {
            return nodes;
        }

        std::string toString() const {
            std::stringstream ss;
            ss << root.symbol() << " | ";
            for (size_t i = 0; i < NSpecies - 1; i++) {
                ss << (i > 0 ? "," : "") << nodes[i].symbol();
            }
            return ss.str();
        }

    private:
        Species root;
        std::array<Species, NSpecies - 1> nodes;
    };
}

#endif