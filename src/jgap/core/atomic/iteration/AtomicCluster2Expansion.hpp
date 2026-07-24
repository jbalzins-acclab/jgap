#ifndef JGAP_ATOMICCLUSTER2EXPANSION_HPP
#define JGAP_ATOMICCLUSTER2EXPANSION_HPP

#include <vector>

#include "Cluster2ExpansionResult.hpp"
#include "jgap/core/CalculationType.hpp"
#include "jgap/core/atomic/geometry/Cluster.hpp"
#include "jgap/core/atomic/neighbours/NeighbourLists.hpp"
#include "jgap/core/atomic/species/composition/Species2Atomic.hpp"
#include "jgap/utils/Utils.hpp"

namespace jgap {

    class AtomicCluster2Expansion {
    public:
        static constexpr size_t ClusterPermutationsAvailable = factorial(1u);

        AtomicCluster2Expansion(const Species2Atomic& species_set) : species_set(species_set) {}

        Cluster2ExpansionResult find(size_t atom_index, const NeighbourLists& neighbour_list,
                                     CalculationType calc_type) const {
            Cluster2ExpansionResult result(calc_type);

            const auto& node_species = species_set.node;
            auto atom_neighbours = neighbour_list.neighbours_per_atom[atom_index].find(node_species);

            if (atom_neighbours != neighbour_list.neighbours_per_atom[atom_index].end()) {
                result.reserve(atom_neighbours->second.size());
                for (auto neigh_data: atom_neighbours->second) {
                    std::array neigh_array{neigh_data};
                    result.add(atom_index, neigh_array);
                }
            }

            return result;
        }

    private:
        Species2Atomic species_set;
    };
}

#endif
