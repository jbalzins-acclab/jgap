#include "TabGapPotential.hpp"

#include "jgap/core/splines/CubicBSpline.hpp"
#include "jgap/core/splines/HermiteCubicSpline.hpp"

namespace jgap {
    TabGapPotential::TabGapPotential(TabulationData energy_tables) {
        isolated_atom_energies = energy_tables.isolated_energies;

        for (auto& [species_pair, grid]: energy_tables.two_body_grids.value_grids) {
            two_body_components.emplace(
                species_pair, TwoBodyTGComponent(
                                  species_pair, energy_tables.eam_grids_vec.empty() ? CubicBSpline::fit(grid)
                                                                                    : HermiteCubicSpline(grid)
                              )
            );
        }

        for (auto& [species_triplet, grid]: energy_tables.three_body_grids.value_grids) {
            three_body_components.emplace(
                species_triplet, ThreeBodyTGComponent(species_triplet, CubicBSpline3D::fit(grid))
            );
        }

        for (auto& eam_grids: energy_tables.eam_grids_vec) {
            EamTGComponent eam_comp(eam_grids);
            Species central = eam_comp.getCentralSpecies();
            eam_components.emplace(central, std::move(eam_comp));
        }
    }

    TabGapPotential::TabGapPotential(
        std::map<Species, Real> isolated_atom_energies,
        std::map<Species2Sorted, TwoBodyTGComponent> two_body_components,
        std::map<Species3AtomicSorted, ThreeBodyTGComponent> three_body_components,
        std::multimap<Species, EamTGComponent> eam_components
    ) :
        isolated_atom_energies(std::move(isolated_atom_energies)),
        two_body_components(std::move(two_body_components)),
        three_body_components(std::move(three_body_components)),
        eam_components(std::move(eam_components)) {}

    AtomicQuantity TabGapPotential::calculateEnergy(const Atoms& atoms) const {
        AtomicQuantity result(atoms.nAtoms());

        static std::set<Species> warned_species;
        for (const Species& species: atoms.getSpecies()) {
            if (isolated_atom_energies.contains(species)) {
                result.value += isolated_atom_energies.at(species);
            } else if (!warned_species.contains(species)) {
                warned_species.insert(species);
                JGAP_LOG_WARN("Isolated atom energy for species {} not found in a tabGAP potential", species.symbol());
            }
        }

        const NeighbourLists nl{atoms, getCutoffs().maxOverall()};
        for (const auto& [pair, comp]: two_body_components) {
            result += comp.energy(nl);
        }
        for (const auto& [triplet, comp]: three_body_components) {
            result += comp.energy(nl);
        }
        for (const auto& [central, comp]: eam_components) {
            result += comp.energy(nl);
        }

        return result;
    }

    Cutoffs TabGapPotential::getCutoffs() const {
        Cutoffs result{};

        for (const auto& [pair, comp]: two_body_components) {
            result += comp.getCutoffs();
        }
        for (const auto& [triplet, comp]: three_body_components) {
            result += comp.getCutoffs();
        }
        for (const auto& [central, comp]: eam_components) {
            result += comp.getCutoffs();
        }

        return result;
    }

    void TabGapPotential::fillTables(TabulationData& table) const {
        for (auto& [species, energy]: isolated_atom_energies) {
            table.isolated_energies[species] += energy;
        }

        for (const auto& [pair, comp]: two_body_components) {
            comp.tabulate(table);
        }
        for (const auto& [triplet, comp]: three_body_components) {
            comp.tabulate(table);
        }
        for (const auto& [central, comp]: eam_components) {
            comp.tabulate(table);
        }
    }
}
