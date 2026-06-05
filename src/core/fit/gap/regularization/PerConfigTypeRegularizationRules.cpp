#include "PerConfigTypeRegularizationRules.hpp"
#include "core/atomic/Atoms.hpp"
#include <cmath>

namespace jgap {
    PerConfigTypeRegularizationRules::PerConfigTypeRegularizationRules(Real energy_sigma_per_atom,
                                                                       Real force_component_sigma,
                                                                       Real virials_iso_sigma_per_atom,
                                                                       Real virials_aniso_sigmas_per_atom,
                                                                       const std::map<std::string, Real> &exact,
                                                                       const std::map<std::string, Real> &contains)
        : default_energy_per_atom(energy_sigma_per_atom),
          exact(exact),
          contains(contains) {

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

    void PerConfigTypeRegularizationRules::fillSigmas(Regularization& sigmas, const Atoms& atoms) const {
        Real multiplier = 1.0;
        const std::string ct = atoms.getConfigType().value_or("default");

        if (exact.contains(ct)) {
            multiplier = exact.at(ct);
        } else {
            std::string longest_match_key;
            for (const auto& [key, value] : contains) {
                if (ct.find(key) != std::string::npos) {
                    if (key.length() > longest_match_key.length()) {
                        longest_match_key = key;
                    }
                }
            }

            if (!longest_match_key.empty()) {
                multiplier = contains.at(longest_match_key);
            }
        }

        sigmas.energy = default_energy_per_atom * multiplier * pow(atoms.nAtoms(), 0.5);
        sigmas.virials = default_virials_per_atom * multiplier * pow(atoms.nAtoms(), 0.5);

        sigmas.forces = std::vector<Vector3>(atoms.nAtoms());
        for (int i = 0; i < atoms.nAtoms(); i++) {
            (*sigmas.forces)[i] = default_force * multiplier * pow(atoms.nAtoms(), 0.5);
        }
    }
}
