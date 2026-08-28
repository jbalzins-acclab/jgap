#include "TwoBodyTGComponent.hpp"
#include "jgap/core/atomic/iteration/Cluster2Expansion.hpp"

namespace jgap {
    TwoBodyTGComponent::TwoBodyTGComponent(Species2Sorted species_pair, ValuePtr<Spline<1>> spline) :
        species_pair(species_pair), spline(std::move(spline)) {}

    AtomicQuantity TwoBodyTGComponent::energy(const NeighbourLists& nl) const {
        AtomicQuantity result(nl.nAtoms());

        Cluster2Expansion expansion(species_pair);
        expansion.forEach(nl, [&](const Cluster2& cluster) {
            auto [E_pair, dE_dr_pair] = spline->interpolate({cluster.separation01.magnitude});

            Real E_cluster = 0.5 * E_pair;
            result.value += E_cluster;

            Real dE_dr = 0.5 * dE_dr_pair[0];
            utils::accumulatePairDistanceDerivatives(
                result.forces[cluster.idx0], result.forces[cluster.idx1], result.virials, dE_dr, cluster.separation01
            );
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
