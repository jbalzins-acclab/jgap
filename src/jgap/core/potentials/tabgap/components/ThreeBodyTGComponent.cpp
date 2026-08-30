#include "ThreeBodyTGComponent.hpp"
#include "jgap/utils/Utils.hpp"

namespace jgap {
    ThreeBodyTGComponent::ThreeBodyTGComponent(const Species3AtomicSorted& species, const CubicBSpline3D& spline) :
        species(species), spline(spline), expansion(species) {}

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
