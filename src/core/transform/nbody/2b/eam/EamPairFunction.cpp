#include "EamPairFunction.hpp"

#include "core/transform/manybody/TwoBodySum.hpp"

namespace jgap {
    std::map<Species, ValuePtr<NBodyAggregator<1>>> EamPairFunction::createAggregators(
        const ValuePtr<TwoBodyTransformation<1>>& base_pf, const std::vector<Atoms>& training_data, EamMode mode) {
        if (auto cast = dynamic_cast<const EamPairFunction*>(base_pf.get()); cast == nullptr) {
            JGAP_LOG_AND_THROW("Non EAM pair function cluster transformation provided");
        }

        std::map<Species, ValuePtr<NBodyAggregator<1>>> aggregators;

        std::set<Species> all_species;
        for (const auto& atoms: training_data) {
            for (const auto& species: atoms.getSpecies()) {
                all_species.insert(species);
            }
        }

        for (const auto& central_species: all_species) {
            auto aggregator = TwoBodySum<1>(central_species);
            auto Z_center_opt = central_species.atomicNumber();
            if (!Z_center_opt && mode != EamMode::Blind && mode != EamMode::EAM) {
                JGAP_LOG_AND_THROW("Central species of unknown atomic number - incompatible with the EAM mode")
            }
            Real Z_center = static_cast<Real>(Z_center_opt.value_or(0));

            for (const auto& contributor_species: all_species) {
                auto pf_clone = base_pf;
                auto& eam_pf_clone = dynamic_cast<EamPairFunction&>(*pf_clone);

                auto Z_contrib_opt = contributor_species.atomicNumber();
                if (!Z_contrib_opt && mode != EamMode::Blind) {
                    JGAP_LOG_AND_THROW("Contributor species of unknown atomic number - incompatible with the EAM mode")
                }

                Real prefactor = 1.0;
                Real Z_contrib = static_cast<Real>(Z_contrib_opt.value_or(0));
                if (mode == EamMode::FSsym) {
                    prefactor = std::sqrt(Z_contrib * Z_center) / 40.0;
                } else if (mode == EamMode::FSgen) {
                    prefactor = std::pow(Z_center, 0.1) * std::sqrt(Z_contrib) / 10.0;
                } else if (mode == EamMode::EAM) {
                    prefactor = std::sqrt(Z_contrib) / 10.0;
                }

                eam_pf_clone.setPrefactor(prefactor);
                aggregator.extend({central_species, contributor_species}, std::move(pf_clone));
            }
            aggregators.emplace(central_species, std::move(aggregator));
        }

        return aggregators;
    }
}
