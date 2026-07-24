#include "ThreeBodyTGComponent.hpp"

namespace jgap {
    ThreeBodyTGComponent::ThreeBodyTGComponent(const Species3AtomicSorted& species, const CubicBSpline3D& spline) :
        species(species), spline(spline), expansion(species) {}

    AtomicQuantity ThreeBodyTGComponent::energy(const NeighbourLists& nl) const {
        AtomicQuantity result(nl.nAtoms());

        auto expansion_result = expansion.findAll(nl, CalculationType::WithGradients);
        if (expansion_result.clusters.empty()) return result;

        assert(expansion_result.derivatives.has_value());

        for (const auto& [cluster, cluster_derivs]:
             std::views::zip(expansion_result.clusters, *expansion_result.derivatives)) {
            auto [q, dq_dr] = transformation.evaluateAndDifferentiate(cluster);
            auto [E, dE_dq] = spline.interpolate(q);

            result.value += E;

            for (size_t i = 0; i < 3; i++) {
                for (size_t j = i + 1; j < 3; j++) {
                    const size_t idx_ij = flattenedIndex(i, j);
                    const Vector3& dir_ij = cluster_derivs.val[idx_ij].direction;
                    const Virials& vir_ij = cluster_derivs.val[idx_ij].virials;

                    Real dE_drij{};

                    for (size_t dim = 0; dim < 3; dim++) {
                        dE_drij += dE_dq[dim] * dq_dr[idx_ij][dim];
                    }

                    result.virials += vir_ij * dE_drij;

                    Vector3 f_ji = dir_ij * dE_drij;
                    result.forces[cluster.atom_indexes[i]] += f_ji;
                    result.forces[cluster.atom_indexes[j]] -= f_ji;
                }
            }
        }

        return result;
    }

    Cutoffs ThreeBodyTGComponent::getCutoffs() const {
        return Cutoffs{{3, std::max(spline.getCutoff()[0], spline.getCutoff()[1])}};
    }

    void ThreeBodyTGComponent::tabulate(TabulationData& tables) const {
        auto& table = tables.three_body_grids.getValueGrid(species);

        for (const auto cell: table) {
            cell.value += spline.interpolate(cell.pos).value;
        }
    }

    std::set<Species> ThreeBodyTGComponent::getAllSpecies() const {
        return {species.root, species.nodes[0], species.nodes[1]};
    }
}
