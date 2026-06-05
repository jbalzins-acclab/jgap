#include "SimpleRegularizationRules.hpp"

namespace jgap {

    SimpleRegularizationRules::SimpleRegularizationRules(Real energy_sigma_per_atom,
                                                         Real force_component_sigma,
                                                         Real virials_iso_sigma_per_atom,
                                                         Real virials_aniso_sigmas_per_atom,
                                                         Real liquid_multiplier,
                                                         Real short_range_multiplier)
        : default_energy_per_atom(energy_sigma_per_atom),
          liquid_multiplier(liquid_multiplier),
          short_range_multiplier(short_range_multiplier) {

        default_force = Vector3{force_component_sigma, force_component_sigma, force_component_sigma};

        default_virials_per_atom = Virials{
            .xx = virials_iso_sigma_per_atom,
            .xy = virials_aniso_sigmas_per_atom,
            .xz = virials_aniso_sigmas_per_atom,

            .yy = virials_iso_sigma_per_atom,
            .yz = virials_aniso_sigmas_per_atom,

            .zz = virials_iso_sigma_per_atom,
        };
    }

    void SimpleRegularizationRules::fillSigmas(Regularization &sigmas,
                                               const Atoms &atoms) const {

        Real multiplier = 1.0;
        const std::string ct = atoms.getConfigType().value_or("default");

        if (ct == "isolated_atom") {
            multiplier = 0.001;
        }

        if (ct.contains("liquid") || ct.contains("melt")) {
            multiplier = liquid_multiplier;
        }

        if (ct.contains("short") || ct.contains("traj") || ct.contains("low_volume")
            || ct.contains("dimer") || ct.contains("trimer")) {
            multiplier = short_range_multiplier;
        }

        sigmas.energy = default_energy_per_atom * multiplier * pow(atoms.nAtoms(), 0.5);
        sigmas.virials = default_virials_per_atom * multiplier * pow(atoms.nAtoms(), 0.5);

        sigmas.forces = std::vector<Vector3>(atoms.nAtoms());
        for (int i = 0; i < atoms.nAtoms(); i++) {
            (*sigmas.forces)[i] = default_force * multiplier * pow(atoms.nAtoms(), 0.5);
        }
    }
}
