#ifndef JGAP_EAMPAIRFUNCTION_HPP
#define JGAP_EAMPAIRFUNCTION_HPP

#include <map>
#include <set>
#include "core/transform/ClusterTransformation.hpp"
#include "core/transform/aggregated/TransformationAggregatorImpl.hpp"
#include "core/potentials/gap/component/ManyBodyGapComponent.hpp"
#include "core/atomic/Atoms.hpp"
#include "core/atomic/neighbours/NeighbourList.hpp"

namespace jgap {
    enum class EamMode {
        FSsym,
        FSgen,
        EAM,
        Blind
    };

    class EamPairFunction : public ClusterTransformation<1, 2> {
    public:
        static std::map<Species, ValuePtr<TransformationAggregator<1>>> createAggregators(
            const ValuePtr<ClusterTransformation<1, 2>>& base_pf,
            const std::vector<Atoms>& training_data,
            EamMode mode) {

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
                if (!Z_center_opt) continue;
                Real Z_center = static_cast<Real>(*Z_center_opt);

                for (const auto& contributor_species : all_species) {
                    auto pf_clone = base_pf;
                    auto& eam_pf_clone = dynamic_cast<EamPairFunction&>(*pf_clone);

                    Real prefactor = 1.0;
                    Real Z_contrib = static_cast<Real>(contributor_species.atomicNumber().value_or(0));
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

        template<CKernelOfDim<1> TKernel>
        static std::vector<ManyBodyGapComponent<1, TKernel>> createComponents(
            const ValuePtr<ClusterTransformation<1, 2>>& base_pf,
            const TKernel& kernel,
            const Sparsifier<1>& sparsifier,
            const std::vector<Atoms>& training_data,
            EamMode mode,
            const std::vector<Real>& optional_coeffs = {}) {

            auto aggregators = createAggregators(base_pf, training_data, mode);
            std::vector<ManyBodyGapComponent<1, TKernel>> components;

            for (auto& [central_species, aggregator] : aggregators) {
                components.emplace_back(std::move(aggregator), kernel, sparsifier, training_data, optional_coeffs);
            }

            return components;
        }

        Cutoffs getCutoffs() const override {
            return Cutoffs{ {2, cutoff} };
        }

        void setPrefactor(Real p) {
            prefactor = p;
        }

        Real getPrefactor() const {
            return prefactor;
        }

        Real getCutoff() const {
            return cutoff;
        }

    protected:
        EamPairFunction(Real cutoff = 0.0, Real prefactor = 1.0)
            : cutoff(cutoff), prefactor(prefactor) {}

        Real cutoff;
        Real prefactor;
    };
}

#endif