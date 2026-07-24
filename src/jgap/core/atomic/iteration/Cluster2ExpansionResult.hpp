#ifndef JGAP_CLUSTER2EXPANSIONRESULT_HPP
#define JGAP_CLUSTER2EXPANSIONRESULT_HPP

#include <optional>
#include <vector>
#include "jgap/core/CalculationType.hpp"
#include "jgap/core/atomic/geometry/Cluster2.hpp"
#include "jgap/core/atomic/geometry/Cluster2Derivatives.hpp"
#include "jgap/core/atomic/neighbours/NeighbourData.hpp"

namespace jgap {

    struct Cluster2ExpansionResult {
        std::vector<Cluster2> clusters;
        std::optional<std::vector<Cluster2Derivatives>> derivatives;

        Cluster2ExpansionResult() = default;

        explicit Cluster2ExpansionResult(CalculationType calc_type) {
            if (calc_type == WithGradients) {
                derivatives.emplace();
            }
        }

        void reserve(size_t n) {
            clusters.reserve(n);
            if (derivatives.has_value()) {
                derivatives->reserve(n);
            }
        }

        void add(size_t index0, const std::array<NeighbourData, 1> &atom_neigh) {
            auto &cluster = clusters.emplace_back();
            cluster.atom_indexes[0] = index0;
            cluster.atom_indexes[1] = atom_neigh[0].neighbour_index;
            cluster.r01 = atom_neigh[0].separation.magnitude;

            if (derivatives.has_value()) {
                auto &deriv = derivatives->emplace_back();
                deriv.dr01 = atom_neigh[0].separation.derivatives;
            }
        }
    };

} // namespace jgap

#endif // JGAP_CLUSTER2EXPANSIONRESULT_HPP
