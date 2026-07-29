#include "TwoBodyTGComponent.hpp"
#include "jgap/core/atomic/iteration/Cluster2Expansion.hpp"

namespace jgap {
    TwoBodyTGComponent::TwoBodyTGComponent(Species2Sorted species_pair, ValuePtr<Spline<1>> spline) :
        species_pair(species_pair), spline(std::move(spline)), expansion(species_pair) {}

    AtomicQuantity TwoBodyTGComponent::energy(const NeighbourLists& nl) const {
        AtomicQuantity result(nl.nAtoms());

        expansion.forEach(nl, [&](const Cluster2& cluster) {
            auto [val, deriv] = spline->interpolate({cluster.separation01.magnitude});

            result.value += val;

            result.virials += cluster.separation01.virials() * deriv[0];

            Vector3 f10 = cluster.separation01.direction * deriv[0];
            result.forces[cluster.idx0] += f10;
            result.forces[cluster.idx1] -= f10;
        });

        return result;
    }

    Cutoffs TwoBodyTGComponent::getCutoffs() const { return Cutoffs{{2u, spline->getCutoff()[0]}}; }

    void TwoBodyTGComponent::tabulate(TabulationData& tables) const {
        auto& table = tables.two_body_grids.getValueGrid(species_pair);

        for (const auto cell: table) {
            cell.value += spline->interpolate(cell.pos).value;
        }
    }
}
