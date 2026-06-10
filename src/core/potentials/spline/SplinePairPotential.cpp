#include "SplinePairPotential.hpp"

#include <vector>
#include <algorithm>

#include "core/atomic/Atoms.hpp"
#include "core/atomic/neighbours/NeighbourList.hpp"
#include "io/log/CurrentLogger.hpp"
#include "utils/Utils.hpp"

namespace jgap {

    void SplinePairPotential::extend(Species species1, Species species2,
                                     const std::vector<Real> &r,
                                     const std::vector<Real> &energies) {

        if (per_species_interpolators.contains(SpeciesSet<2, Symmetric>{species1, species2})) {
            JGAP_LOG_AND_THROW("Trying to re-write interpolators for {}-{} pair", species1.symbol(), species2.symbol());
        }

        per_species_interpolators.emplace(
            SpeciesSet<2, Symmetric>{species1, species2}, NaturalCubicSpline(r, energies)
            );
    }

    AtomicQuantity SplinePairPotential::calculateEnergy(const Atoms &atoms) {

        AtomicQuantity result(atoms.nAtoms());

        NeighbourList nl(atoms, getCutoffs().maxOverall());

        for (const auto &[species_pair, interpolator]: per_species_interpolators) {

            auto pair_clusters = nl.findAllClusters(species_pair);

            for (auto pair: pair_clusters) {

                const Separation& separation = pair.between(0, 1);

                InterpolationResults<1> spline_val = interpolator.interpolate({separation.magnitude});

                result.value += spline_val.value * 0.5;

                Real dE_dr = spline_val.gradient[0] * 0.5;

                result.virials += separation.virials * dE_dr;

                Vector3 f01 = separation.direction * dE_dr;
                result.forces[pair.atom_indexes[1]] -= f01;
                result.forces[pair.atom_indexes[0]] += f01;
            }
        }

        return result;
    }

    Cutoffs SplinePairPotential::getCutoffs() {
        double cutoff = 0.0;
        for (const auto& interpolator: per_species_interpolators | std::views::values) {
            cutoff = std::max(cutoff, interpolator.getCutoff()[0]);
        }
        return {{2, cutoff}};
    }
}
