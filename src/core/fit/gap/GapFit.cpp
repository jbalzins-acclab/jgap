#include "GapFit.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {
    void GapFit::fit(GapPotential &to_be_fit,
                     const std::vector<Atoms> &training_data,
                     const std::vector<NeighbourList> &neighbour_lists,
                     const std::vector<Regularization> &sigmas,
                     const std::shared_ptr<RegularizationRules> &regularization_rules) {

        auto sigmas_inverse = sigmas;
        if (sigmas_inverse.empty()) {
            sigmas_inverse.resize(training_data.size());
        }
        if (sigmas_inverse.size() != training_data.size()) {
            JGAP_LOG_AND_THROW("Number of training structures {} doesn't match number of regularization params {}",
                training_data.size(), sigmas_inverse.size());
        }

        if (regularization_rules) {
            for (size_t i = 0; i < training_data.size(); i++) {
                regularization_rules->fillSigmas(sigmas_inverse[i], training_data[i], neighbour_lists[i]);

                sigmas_inverse[i].energy = sigmas_inverse[i].energy.transform(
                    [](Real val) -> Real {return 1.0 / val;}
                    );
                sigmas_inverse[i].virials = sigmas_inverse[i].virials.transform([](const Virials &val) -> Virials {
                    return Virials{
                        static_cast<Real>(1.0) / val.xx,
                        static_cast<Real>(1.0)  / val.xy,
                        static_cast<Real>(1.0)  / val.xz,
                        static_cast<Real>(1.0)  / val.yy,
                        static_cast<Real>(1.0)  / val.yz,
                        static_cast<Real>(1.0)  / val.zz
                    };
                });
                sigmas_inverse[i].forces = sigmas_inverse[i].forces.transform(
                    [](const std::vector<Vector3>& val) -> std::vector<Vector3> {

                    std::vector<Vector3> result(val.size());
                    for (size_t j = 0; j < val.size(); j++) {
                        result[j] = Vector3{
                            static_cast<Real>(1.0) / val[j].x,
                            static_cast<Real>(1.0) / val[j].y,
                            static_cast<Real>(1.0) / val[j].z
                        };
                    }

                    return result;
                });
            }
        }

        Real required_cutoff = 0.0;

        for (const auto& component : to_be_fit.components) {
            required_cutoff = std::max(required_cutoff, component->getCutoffs().maxOverall());
        }

        if (to_be_fit.optional_external_potential) {
            required_cutoff = std::max(
                required_cutoff,
                to_be_fit.optional_external_potential->getCutoffs().maxOverall()
                );
        }

        std::vector<NeighbourList> actual_neighbour_lists;

        if (neighbour_lists.empty()) {
            actual_neighbour_lists = NeighbourList::form(training_data, required_cutoff);
        } else {

            if (neighbour_lists.size() != training_data.size()) {
                JGAP_LOG_AND_THROW("Number of training structures {} doesn't match number of neighbour lists {}.",
                    training_data.size(), neighbour_lists.size());
            }

            actual_neighbour_lists = neighbour_lists;
            for (const auto& nl : actual_neighbour_lists) {
                if (nl.cutoff < required_cutoff) {
                    JGAP_LOG_AND_THROW("Provided neighbour list cutoff ({}) is smaller than the required cutoff ({}).",
                        nl.cutoff, required_cutoff);
                }
            }
        }

        std::vector<EnergyData> energies_without_external(training_data.size());

        for (size_t i = 0; i < training_data.size(); ++i) {
            const auto& atoms = training_data[i];

            EnergyData data;

            data.energy = atoms.getEnergy();
            data.virials = atoms.getVirials();
            data.forces = atoms.getForces();

            energies_without_external[i] = data;
        }

        if (to_be_fit.optional_external_potential) {
            for (size_t i = 0; i < training_data.size(); ++i) {
                // Here we pass the original Atoms object to the external potential
                // since it's required by Potential::energy
                AtomicQuantity ext_pred = to_be_fit.optional_external_potential->energy(training_data[i]);

                if (energies_without_external[i].energy.has_value()) {
                    energies_without_external[i].energy.value() -= ext_pred.value;
                }

                if (energies_without_external[i].virials.has_value()) {
                    energies_without_external[i].virials.value() -= ext_pred.virials;
                }

                if (energies_without_external[i].forces.has_value()) {
                    auto& forces = energies_without_external[i].forces.value();
                    for (size_t atom_idx = 0; atom_idx < forces.size(); ++atom_idx) {
                        forces[atom_idx] -= ext_pred.forces[atom_idx];
                    }
                }
            }
        }

        auto fit_coefficients = mainFit(
            to_be_fit.components,
            actual_neighbour_lists,
            energies_without_external,
            sigmas_inverse
            );
        to_be_fit.setCoefficients(fit_coefficients);

        for (Real coeff: fit_coefficients) {
            JGAP_LOG_INFO("Coefficient: {}", coeff);
        }
    }
}
