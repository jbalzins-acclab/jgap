#include "EamTGComponent.hpp"

#include "../../../transform/eam/SplinePairTransformation.hpp"
#include "core/transform/aggregated/TransformationAggregatorImpl.hpp"

namespace jgap {

    EamTGComponent::EamTGComponent(ManyBodyGrids<2> grids_for_element)
        : energy_spline(HermiteCubicSpline{std::move(grids_for_element.value_grid)}),
          spline_density_aggregator{grids_for_element.central_atom_species} {

        for (auto &[species_pair, eam_pf_grid]: grids_for_element.aggregator_grids.value_grids) {
            spline_density_aggregator.extend(
                species_pair,
                SplinePairTransformation(HermiteCubicSpline(std::move(eam_pf_grid)))
            );
        }
    }

    AtomicQuantity EamTGComponent::energy(const NeighbourList &nl) const {
        size_t n = nl.nAtoms();
        AtomicQuantity result(n);

        auto descriptors = spline_density_aggregator.aggregate(nl);
        for (const auto& [atom_idx, atom_descriptor]: descriptors) {
            auto [E, dE_dq] = energy_spline.interpolate(atom_descriptor.value);

            result.value += E;
            result.virials += atom_descriptor.virials[0] * dE_dq[0];

            for (size_t i = 0; i < n; i++) {
                result.forces[i] += atom_descriptor.forces[i][0] * dE_dq[0];
            }
        }

        return result;
    }

    Cutoffs EamTGComponent::getCutoffs() const {
        return spline_density_aggregator.getCutoffs();
    }

    void EamTGComponent::tabulate(TabulationData &tables) const {

        spline_density_aggregator.tabulateNewManyBodyGrid(tables);
        auto& value_grid = tables.eam_grids_vec.back().value_grid;

        for (auto cell: value_grid) {
            cell.value += energy_spline.interpolate(cell.pos).value;
        }
    }

    std::set<Species> EamTGComponent::getAllSpecies() const {
        return spline_density_aggregator.getAllSpecies();
    }
}
