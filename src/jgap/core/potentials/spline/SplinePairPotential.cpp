#include "SplinePairPotential.hpp"

#include <algorithm>
#include <vector>

#include "../../io/log/CurrentLogger.hpp"
#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/core/atomic/iteration/Cluster2Expansion.hpp"
#include "jgap/core/atomic/neighbours/NeighbourLists.hpp"
#include "jgap/utils/Utils.hpp"

namespace jgap {

    void SplinePairPotential::extend(
        Species species1, Species species2, const std::vector<Real>& r, const std::vector<Real>& energies
    ) {
        if (per_species_interpolators.contains(Species2Sorted{species1, species2})) {
            JGAP_LOG_AND_THROW("Trying to re-write interpolators for {}-{} pair", species1.symbol(), species2.symbol());
        }

        per_species_interpolators.emplace(Species2Sorted{species1, species2}, NaturalCubicSpline(r, energies));
    }

    AtomicQuantity SplinePairPotential::calculateEnergy(const Atoms& atoms) const {
        AtomicQuantity result(atoms.nAtoms());

        NeighbourLists nl(atoms, getCutoffs().maxOverall());

        for (const auto& [species_pair, interpolator]: per_species_interpolators) {
            Cluster2Expansion expansion(species_pair);
            expansion.forEach(nl, [&](const Cluster2& cluster) {
                auto [E_pair, gradient] = interpolator.interpolate({cluster.separation01.magnitude});

                Real E_cluster = 0.5 * E_pair;
                result.value += E_cluster;

                Real dE_dr = 0.5 * gradient[0];

                utils::accumulatePairDistanceDerivatives(
                    result.forces[cluster.idx0],
                    result.forces[cluster.idx1],
                    result.virials,
                    dE_dr,
                    cluster.separation01
                );
            });
        }

        return result;
    }

    void SplinePairPotential::fillTables(TabulationData& tables) const {
        for (const auto& [species_pair, interpolator]: per_species_interpolators) {
            auto& table = tables.two_body_grids.getValueGrid(species_pair);
            for (const auto cell: table) {
                cell.value += interpolator.interpolate(cell.pos).value;
            }
        }
    }

    Cutoffs SplinePairPotential::getCutoffs() const {
        Real cutoff = 0.0;
        for (const auto& interpolator: per_species_interpolators | std::views::values) {
            cutoff = std::max(cutoff, interpolator.getCutoff()[0]);
        }
        return {{2, cutoff}};
    }
}
