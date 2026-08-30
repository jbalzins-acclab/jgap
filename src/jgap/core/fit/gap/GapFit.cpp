#include "GapFit.hpp"
#include "../../io/log/CurrentLogger.hpp"
#include "jgap/core/UnseqFor.hpp"

namespace jgap {
    void GapFit::fit(
        GapPotential& to_be_fit, const std::vector<Atoms>& training_data, const std::vector<Regularization>& sigmas
    ) {
        if (sigmas.size() != training_data.size()) {
            JGAP_LOG_AND_THROW(
                "Number of training structures ({}) doesn't match number of regularization sigmas ({})",
                training_data.size(),
                sigmas.size()
            );
        }

        auto sigmas_inverse = sigmas;

        for (size_t i = 0; i < training_data.size(); i++) {
            const auto& atoms = training_data[i];
            const auto& reg = sigmas[i];

            // Validation checks: if data is present, regularization must be specified and positive
            if (atoms.getEnergy().has_value()) {
                if (!reg.energy.has_value() || *reg.energy <= 0.0_r) {
                    JGAP_LOG_AND_THROW(
                        "Structure {} has energy data, but no valid positive energy regularization sigma was specified",
                        i
                    );
                }
            }

            if (atoms.getForces().has_value()) {
                if (!reg.forces.has_value() || reg.forces->size() != atoms.nAtoms()) {
                    JGAP_LOG_AND_THROW(
                        "Structure {} has forces data (size {}), but force regularization sigmas are missing or size "
                        "mismatch",
                        i,
                        atoms.nAtoms()
                    );
                }
                for (size_t a = 0; a < atoms.nAtoms(); ++a) {
                    const auto& f_sig = (*reg.forces)[a];
                    if (f_sig.x <= 0.0_r || f_sig.y <= 0.0_r || f_sig.z <= 0.0_r) {
                        JGAP_LOG_AND_THROW("Structure {} has non-positive force regularization sigma at atom {}", i, a);
                    }
                }
            }

            if (atoms.getVirials().has_value()) {
                if (!reg.virials.has_value()) {
                    JGAP_LOG_AND_THROW(
                        "Structure {} has virials data, but no virial regularization sigmas were specified", i
                    );
                }
                const auto& v_sig = *reg.virials;
                if (v_sig.xx <= 0.0_r || v_sig.yy <= 0.0_r || v_sig.zz <= 0.0_r || v_sig.xy <= 0.0_r
                    || v_sig.xz <= 0.0_r || v_sig.yz <= 0.0_r) {
                    JGAP_LOG_AND_THROW("Structure {} has non-positive virial regularization sigma", i);
                }
            }

            sigmas_inverse[i].energy = sigmas_inverse[i].energy.transform([](Real val) -> Real { return 1.0_r / val; });
            sigmas_inverse[i].virials = sigmas_inverse[i].virials.transform([](const Virials& val) -> Virials {
                return Virials{
                    1.0_r / val.xx, 1.0_r / val.xy, 1.0_r / val.xz, 1.0_r / val.yy, 1.0_r / val.yz, 1.0_r / val.zz
                };
            });
            sigmas_inverse[i].forces =
                sigmas_inverse[i].forces.transform([](const std::vector<Vector3>& val) -> std::vector<Vector3> {
                    std::vector<Vector3> result(val.size());
                    for (size_t j = 0; j < val.size(); j++) {
                        result[j] = Vector3{1.0_r / val[j].x, 1.0_r / val[j].y, 1.0_r / val[j].z};
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
