#include "GapFit.hpp"
#include "jgap/io/log/CurrentLogger.hpp"

#include "jgap/core/UnseqFor.hpp"

namespace jgap {
    void GapFit::fit(GapPotential& to_be_fit, const std::vector<Atoms>& training_data,
                     const RegularizationRules& regularization_rules, const std::vector<Regularization>& sigmas) {

        auto sigmas_inverse = sigmas;
        if (sigmas_inverse.empty()) {
            sigmas_inverse.resize(training_data.size());
        }
        if (sigmas_inverse.size() != training_data.size()) {
            JGAP_LOG_AND_THROW("Number of training structures {} doesn't match number of regularization params {}",
                               training_data.size(), sigmas_inverse.size());
        }

        for (size_t i = 0; i < training_data.size(); i++) {
            regularization_rules.fillSigmas(sigmas_inverse[i], training_data[i]);

            sigmas_inverse[i].energy = sigmas_inverse[i].energy.transform([](Real val) -> Real { return 1.0 / val; });
            sigmas_inverse[i].virials = sigmas_inverse[i].virials.transform([](const Virials& val) -> Virials {
                return Virials{1.0 / val.xx, 1.0 / val.xy, 1.0 / val.xz, 1.0 / val.yy, 1.0 / val.yz, 1.0 / val.zz};
            });
            sigmas_inverse[i].forces =
                sigmas_inverse[i].forces.transform([](const std::vector<Vector3>& val) -> std::vector<Vector3> {
                    std::vector<Vector3> result(val.size());
                    for (size_t j = 0; j < val.size(); j++) {
                        result[j] = Vector3{1.0 / val[j].x, 1.0 / val[j].y, 1.0 / val[j].z};
                    }

                    return result;
                });
        }

        std::vector<EnergyData> energies_without_external(training_data.size());

        unseqForIndex(0, training_data.size(), [&](size_t i) {
            const auto& atoms = training_data[i];
            auto& energy_without_external = energies_without_external[i];

            energy_without_external = {atoms.getEnergy(), atoms.getVirials(), atoms.getForces()};

            if (to_be_fit.optional_external_potential) {
                AtomicQuantity ext_pred = to_be_fit.optional_external_potential->calculateEnergy(atoms);

                if (energy_without_external.energy.has_value()) {
                    energy_without_external.energy.value() -= ext_pred.value;
                }

                if (energy_without_external.virials.has_value()) {
                    energy_without_external.virials.value() -= ext_pred.virials;
                }

                if (energy_without_external.forces.has_value()) {
                    auto& forces = energy_without_external.forces.value();
                    for (size_t atom_idx = 0; atom_idx < forces.size(); ++atom_idx) {
                        forces[atom_idx] -= ext_pred.forces[atom_idx];
                    }
                }
            }
        });

        auto fit_coefficients =
            findCoefficients(to_be_fit.components, training_data, energies_without_external, sigmas_inverse);
        to_be_fit.setCoefficients(fit_coefficients);
    }
}
