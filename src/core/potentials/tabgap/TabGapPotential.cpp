#include "TabGapPotential.hpp"

#include "components/ThreeBodyTGComponent.hpp"
#include "components/TwoBodyTGComponent.hpp"
#include "core/transform/aggregated/TransformationAggregatorImpl.hpp"

namespace jgap {
    TabGapPotential::TabGapPotential(TabulationData energy_tables) {

        isolated_atom_energies = energy_tables.isolated_energies;

        for (auto& [species_pair, grid]: energy_tables.two_body_grids.value_grids) {
            ++n_2b_components;

            components.emplace_back(TwoBodyTGComponent(
                species_pair,
                energy_tables.eam_grids_vec.empty() ?
                    CubicBSpline::fit(grid) :
                    HermiteCubicSpline(grid)
                ));
        }

        for (auto& [species_triplet, grid]: energy_tables.three_body_grids.value_grids) {
            ++n_3b_components;

            components.emplace_back(ThreeBodyTGComponent(
                species_triplet,
                CubicBSpline3D::fit(grid)
                ));
        }

        for (auto& eam_grids: energy_tables.eam_grids_vec) {
            ++n_eam_components;

            components.emplace_back(EamTGComponent(eam_grids));
        }
    }

    TabGapPotential::TabGapPotential(const std::map<Species, Real> &isolated_atom_energies,
                                     const std::vector<ValuePtr<TabGapComponent> > &components)
        : isolated_atom_energies(isolated_atom_energies),
          components(components) {
        recomputeComponentCounts();
    }

    void TabGapPotential::recomputeComponentCounts() {
        n_2b_components = 0;
        n_3b_components = 0;
        n_eam_components = 0;
        for (const auto& component: components) {
            if (dynamic_cast<const TwoBodyTGComponent*>(component.get())) {
                ++n_2b_components;
            } else if (dynamic_cast<const ThreeBodyTGComponent*>(component.get())) {
                ++n_3b_components;
            } else if (dynamic_cast<const EamTGComponent*>(component.get())) {
                ++n_eam_components;
            }
        }
    }

    AtomicQuantity TabGapPotential::calculateEnergy(const Atoms &atoms) const {

        AtomicQuantity result(atoms.nAtoms());

        for (const Species& species: atoms.lookupSpecies()) {
            if (isolated_atom_energies.contains(species)) {
                result.value += isolated_atom_energies.at(species);
            } else {
                JGAP_LOG_WARN("Unknown species {} for a tabGAP potential;"
                              "comment out this line if that's expected", species.symbol());
            }
        }

        const NeighbourList nl{atoms, getCutoffs().maxOverall()};
        for (const auto& component: components) {
            result += component->energy(nl);
        }

        return result;
    }

    Cutoffs TabGapPotential::getCutoffs() const {
        Cutoffs result{};

        for (const auto& component: components) {
            result += component->getCutoffs();
        }

        return result;
    }

    void TabGapPotential::fillTables(TabulationData &table) const {

        for (auto& [species, energy]: isolated_atom_energies) {
            table.isolated_energies[species] += energy;
        }

        for (const auto& component: components) {
            component->tabulate(table);
        }
    }
}
