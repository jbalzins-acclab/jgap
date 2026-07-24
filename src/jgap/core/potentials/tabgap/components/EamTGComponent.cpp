#include "EamTGComponent.hpp"
#include "jgap/core/transform/nbody/2b/eam/SplinePairTransformation.hpp"


namespace jgap {

    EamTGComponent::EamTGComponent(ManyBodyGrids2<1, 1> grids_for_element) :
        energy_spline(HermiteCubicSpline{grids_for_element.value_grid}),
        spline_density_aggregator{grids_for_element.central_atom_species} {
        for (auto& [species_pair, eam_pf_grid]: grids_for_element.aggregator_grids.value_grids) {
            spline_density_aggregator.extend(species_pair,
                                             SplinePairTransformation(HermiteCubicSpline(std::move(eam_pf_grid))));
        }
    }

    AtomicQuantity EamTGComponent::energy(const NeighbourLists& nl) const {
        size_t n = nl.nAtoms();
        AtomicQuantity result(n);

        auto descriptors = spline_density_aggregator.aggregate(nl);
        for (size_t i = 0; i < descriptors.values.size(); i++) {
            auto [E, dE_dq] = energy_spline.interpolate(descriptors.values[i]);

            result.value += E;
            result.virials += descriptors.virials[i][0] * dE_dq[0];

            for (size_t j = 0; j < n; j++) {
                result.forces[j] += descriptors.forces[i * n + j][0] * dE_dq[0];
            }
        }

        return result;
    }

    Cutoffs EamTGComponent::getCutoffs() const { return spline_density_aggregator.getCutoffs(); }

    void EamTGComponent::tabulate(TabulationData& tables) const {
        spline_density_aggregator.tabulateNewManyBodyGrid(tables);
        auto& value_grid = tables.eam_grids_vec.back().value_grid;

        for (auto cell: value_grid) {
            cell.value += energy_spline.interpolate(cell.pos).value;
        }
    }

    std::set<Species> EamTGComponent::getAllSpecies() const { return spline_density_aggregator.getAllSpecies(); }
}
