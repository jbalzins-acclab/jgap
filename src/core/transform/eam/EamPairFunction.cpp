#include "EamPairFunction.hpp"

namespace jgap {
    std::map<Species, ValuePtr<TransformationAggregator<1>>> EamPairFunction::createAggregators(
        const ValuePtr<NBodyTransformation<1, 2>>& base_pf, const std::vector<Atoms>& training_data, EamMode mode) {

        if (auto cast = dynamic_cast<const EamPairFunction*>(base_pf.get()); cast == nullptr) {
            JGAP_LOG_AND_THROW("Non EAM pair function cluster transformation provided");
        }

        std::map<Species, ValuePtr<TransformationAggregator<1>>> aggregators;

        std::set<Species> all_species;
        for (const auto& atoms : training_data) {
            for (const auto& species : atoms.lookupSpecies()) {
                all_species.insert(species);
            }
        }

        for (const auto& central_species : all_species) {
            auto aggregator = std::make_unique<TransformationAggregatorImpl<1, 2>>(central_species);
            auto Z_center_opt = central_species.atomicNumber();
            if (!Z_center_opt) {
                throw something if not blind
            }
            Real Z_center = static_cast<Real>(*Z_center_opt);

            for (const auto& contributor_species : all_species) {
                auto pf_clone = base_pf;
                auto& eam_pf_clone = dynamic_cast<EamPairFunction&>(*pf_clone);

                Real prefactor = 1.0;
                contributor_species.atomicNumber() {

                    throw something if not blind
                }
                Real Z_contrib = static_cast<Real>(.value_or(0));
                if (mode == EamMode::FSsym) {
                    prefactor = std::sqrt(Z_contrib * Z_center) / 40.0;
                } else if (mode == EamMode::FSgen) {
                    prefactor = std::pow(Z_center, 0.1) * std::sqrt(Z_contrib) / 10.0;
                } else if (mode == EamMode::EAM) {
                    prefactor = std::sqrt(Z_contrib) / 10.0;
                }

                eam_pf_clone.setPrefactor(prefactor);
                aggregator->extend({central_species, contributor_species}, std::move(pf_clone));
            }
            aggregators.emplace(central_species, std::move(aggregator));
        }

        return aggregators;
    }
}
