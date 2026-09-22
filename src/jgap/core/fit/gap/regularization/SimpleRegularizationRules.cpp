#include "SimpleRegularizationRules.hpp"
#include <cmath>
#include "jgap/core/atomic/Atoms.hpp"

namespace jgap {

    SimpleRegularizationRules::SimpleRegularizationRules(
        double energy_sigma_per_atom,
        double force_component_sigma,
        double virials_iso_sigma_per_atom,
        double virials_aniso_sigmas_per_atom,
        double liquid_multiplier,
        double short_range_multiplier
    ) :
        defaults(
            energy_sigma_per_atom, force_component_sigma, virials_iso_sigma_per_atom, virials_aniso_sigmas_per_atom
        ),
        liquid_multiplier(liquid_multiplier),
        short_range_multiplier(short_range_multiplier) {}

    Regularization SimpleRegularizationRules::determine(const Atoms& atoms) const {
        double multiplier = 1.0;
        const std::string ct = atoms.getConfigType().value_or("default");

        if (ct == "isolated_atom") {
            multiplier = 0.001;
        }

        if (ct.contains("liquid") || ct.contains("melt")) {
            multiplier = liquid_multiplier;
        }

        if (ct.contains("short") || ct.contains("traj") || ct.contains("low_volume") || ct.contains("dimer")
            || ct.contains("trimer")) {
            multiplier = short_range_multiplier;
        }

        Regularization sigmas;
        sigmas.energy = defaults.energy * multiplier * std::sqrt(static_cast<double>(atoms.nAtoms()));
        sigmas.virials = defaults.virials * multiplier * std::sqrt(static_cast<double>(atoms.nAtoms()));

        sigmas.forces = std::vector<Vector3>(atoms.nAtoms());
        for (size_t i = 0; i < atoms.nAtoms(); i++) {
            (*sigmas.forces)[i] = defaults.force * multiplier;
        }
        return sigmas;
    }
}
