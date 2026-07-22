#include "SplinePairPotential.hpp"

#include <algorithm>
#include <vector>

#include "core/atomic/Atoms.hpp"
#include "core/atomic/iteration/Cluster2Expansion.hpp"
#include "core/atomic/neighbours/NeighbourLists.hpp"
#include "io/log/CurrentLogger.hpp"
#include "utils/Utils.hpp"

namespace jgap {

    void SplinePairPotential::extend(Species species1, Species species2, const std::vector<Real>& r,
                                     const std::vector<Real>& energies) {
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
            auto expansion_result = expansion.find(nl, CalculationType::WithGradients);

            assert(expansion_result.derivatives.has_value());

            for (const auto& [cluster, cluster_derivs]:
                 std::views::zip(expansion_result.clusters, *expansion_result.derivatives)) {
                const Real& separation_magnitude = cluster.r01;
                const SeparationDerivatives& separation_deriv = cluster_derivs.dr01;

                InterpolationResults<1> spline_val = interpolator.interpolate({separation_magnitude});

                result.value += spline_val.value;

                Real dE_dr = spline_val.gradient[0];

                result.virials += separation_deriv.virials * dE_dr;

                Vector3 f01 = separation_deriv.direction * dE_dr;
                result.forces[cluster.atom_indexes[1]] -= f01;
                result.forces[cluster.atom_indexes[0]] += f01;
            }
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
        double cutoff = 0.0;
        for (const auto& interpolator: per_species_interpolators | std::views::values) {
            cutoff = std::max(cutoff, interpolator.getCutoff()[0]);
        }
        return {{2, cutoff}};
    }
}
