#include "ThreeBodyTGComponent.hpp"
#include "jgap/utils/Utils.hpp"

namespace jgap {
    ThreeBodyTGComponent::ThreeBodyTGComponent(const Species3AtomicSorted& species, const CubicBSpline3D& spline) :
        species(species), spline(spline), expansion(species, ClusterPermutationMode::NoNodePermutation) {
        if (species.nodes[0] == species.nodes[1]) {
            const auto& coeff_grid = spline.getCoefficients();
            if (coeff_grid.sizes[0] != coeff_grid.sizes[1]) {
                JGAP_LOG_AND_THROW(
                    "3-body spline for symmetric species {}-{}-{} must have matching r1 and r2 grid sizes ({} vs {})",
                    species.root.symbol(),
                    species.nodes[0].symbol(),
                    species.nodes[1].symbol(),
                    coeff_grid.sizes[0],
                    coeff_grid.sizes[1]
                );
            }
            if (std::abs(coeff_grid.spacing[0] - coeff_grid.spacing[1]) > 1e-10
                || std::abs(coeff_grid.origin[0] - coeff_grid.origin[1]) > 1e-10) {
                JGAP_LOG_AND_THROW(
                    "3-body spline for symmetric species {}-{}-{} must have matching r1 and r2 grid bounds/spacings",
                    species.root.symbol(),
                    species.nodes[0].symbol(),
                    species.nodes[1].symbol()
                );
            }
            for (size_t i = 0; i < coeff_grid.sizes[0]; ++i) {
                for (size_t j = i + 1; j < coeff_grid.sizes[1]; ++j) {
                    for (size_t k = 0; k < coeff_grid.sizes[2]; ++k) {
                        Real v1 = coeff_grid({i, j, k});
                        Real v2 = coeff_grid({j, i, k});
                        if (std::abs(v1 - v2) > 1e-4) {
                            JGAP_LOG_AND_THROW(
                                "3-body spline for symmetric species {}-{}-{} is not symmetric in swapping distances "
                                "at indices ({}, {}, {}): {} vs {}",
                                species.root.symbol(),
                                species.nodes[0].symbol(),
                                species.nodes[1].symbol(),
                                i,
                                j,
                                k,
                                v1,
                                v2
                            );
                        }
                    }
                }
            }
        }
    }

    AtomicQuantity ThreeBodyTGComponent::energy(const NeighbourLists& nl) const {
        AtomicQuantity result(nl.nAtoms());

        expansion.forEach(nl, [&](const Cluster3& cluster) {
            auto desc = transformation.evaluateAndDifferentiate(cluster);
            auto [E, dE_dq] = spline.interpolate(desc.value);

            result.value += E;

            Vector3 f1{0.0_r, 0.0_r, 0.0_r};
            Vector3 f2{0.0_r, 0.0_r, 0.0_r};
            for (size_t dim = 0; dim < 3; dim++) {
                f1 -= dE_dq[dim] * desc.grad_r1[dim];
                f2 -= dE_dq[dim] * desc.grad_r2[dim];
            }

            result.forces[cluster.idx1] += f1;
            result.forces[cluster.idx2] += f2;
            result.forces[cluster.idx0] -= (f1 + f2);

            result.virials += Virials::dyadic(cluster.separation01.vec(), f1);
            result.virials += Virials::dyadic(cluster.separation02.vec(), f2);
        });

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
