#ifndef JGAP_CLUSTER3EXPANSIONRESULT_HPP
#define JGAP_CLUSTER3EXPANSIONRESULT_HPP

#include <optional>
#include <vector>
#include "jgap/core/CalculationType.hpp"
#include "jgap/core/atomic/geometry/Cluster3.hpp"
#include "jgap/core/atomic/geometry/Cluster3Derivatives.hpp"
#include "jgap/core/atomic/neighbours/NeighbourData.hpp"

namespace jgap {

    struct Cluster3ExpansionResult {
        std::vector<Cluster3> clusters;
        std::optional<std::vector<Cluster3Derivatives>> derivatives;

        Cluster3ExpansionResult() = default;

        explicit Cluster3ExpansionResult(CalculationType calc_type) {
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

        void add(size_t index0, const std::array<NeighbourData, 2>& atom_neigh) {
            auto& cluster = clusters.emplace_back();
            cluster.atom_indexes[0] = index0;
            cluster.atom_indexes[1] = atom_neigh[0].neighbour_index;
            cluster.atom_indexes[2] = atom_neigh[1].neighbour_index;

            cluster.separation_magnitudes[0] = atom_neigh[0].separation.magnitude;
            cluster.separation_magnitudes[1] = atom_neigh[1].separation.magnitude;

            Vector3 r01 = atom_neigh[0].separation.vec();
            Vector3 r02 = atom_neigh[1].separation.vec();

            if (derivatives.has_value()) {
                auto separation12 = Separation(r01, r02);
                cluster.separation_magnitudes[2] = separation12.magnitude;

                auto& deriv = derivatives->emplace_back();
                deriv.val[0] = atom_neigh[0].separation.derivatives;
                deriv.val[1] = atom_neigh[1].separation.derivatives;
                deriv.val[2] = separation12.derivatives;
            } else {
                cluster.separation_magnitudes[2] = (r01 - r02).norm();
            }
        }
    };

}

#endif
