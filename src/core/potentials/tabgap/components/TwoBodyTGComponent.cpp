#include "TwoBodyTGComponent.hpp"

namespace jgap {
    TwoBodyTGComponent::TwoBodyTGComponent(SpeciesSet<2, Symmetric> species_pair, ValuePtr<Spline<1>> spline)
        : species_pair(species_pair), spline(std::move(spline)) {
    }

    AtomicQuantity TwoBodyTGComponent::energy(const NeighbourList &nl) const {
        AtomicQuantity result(nl.nAtoms());

        for (const auto& cluster: nl.findAllClusters<WithDerivatives>(species_pair)) {

            auto [val, deriv] = spline->interpolate({cluster.between(0, 1)});

            result.value += val * SymmetryFactor;

            result.virials += cluster.derivatives[0].virials * deriv[0] * SymmetryFactor;

            Vector3 f10 = cluster.derivatives[0].direction * deriv[0];
            result.forces[cluster.atom_indexes[0]] += f10;
            result.forces[cluster.atom_indexes[1]] -= f10;
        }

        return result;
    }

    Cutoffs TwoBodyTGComponent::getCutoffs() const {
        return Cutoffs{
            {2u, spline->getCutoff()[0]}
        };
    }

    void TwoBodyTGComponent::tabulate(TabulationData &tables) const {
        auto& table = tables.two_body_grids.getValueGrid(species_pair);

        for (const auto cell: table) {
            cell.value += spline->interpolate(cell.pos).value;
        }
    }
}
