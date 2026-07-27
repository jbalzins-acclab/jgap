#include "TwoBodyTGComponent.hpp"
#include "jgap/core/atomic/iteration/Cluster2Expansion.hpp"

namespace jgap {
    TwoBodyTGComponent::TwoBodyTGComponent(Species2Sorted species_pair, ValuePtr<Spline<1>> spline) :
        species_pair(species_pair), spline(std::move(spline)), expansion(species_pair) {}

    AtomicQuantity TwoBodyTGComponent::energy(const NeighbourLists& nl) const {
        AtomicQuantity result(nl.nAtoms());

        auto expansion_result = expansion.find(nl, CalculationType::WithGradients);
        assert(expansion_result.derivatives.has_value());

        for (const auto& [cluster, cluster_derivs]:
             std::views::zip(expansion_result.clusters, *expansion_result.derivatives)) {
            auto [val, deriv] = spline->interpolate({cluster.r01});

            result.value += val;

            result.virials += cluster_derivs.dr01.virials * deriv[0];

            Vector3 f10 = cluster_derivs.dr01.direction * deriv[0];
            result.forces[cluster.atom_indexes[0]] += f10;
            result.forces[cluster.atom_indexes[1]] -= f10;
        }

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
