#ifndef JGAP_GAPCOMPONENTUTILS_HPP
#define JGAP_GAPCOMPONENTUTILS_HPP

#include <map>
#include <optional>
#include <set>
#include <vector>

#include "jgap/core/UnseqFor.hpp"
#include "jgap/core/ValuePtr.hpp"
#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/core/atomic/neighbours/NeighbourLists.hpp"
#include "jgap/core/potentials/gap/component/ManyBodyGapComponent.hpp"
#include "jgap/core/potentials/gap/component/ThreeBodyGapComponent.hpp"
#include "jgap/core/potentials/gap/component/TwoBodyGapComponent.hpp"
#include "jgap/core/transform/manybody/TwoBodySum.hpp"
#include "jgap/core/transform/nbody/2b/eam/EamPairFunction.hpp"
#include "jgap/experimental/transform/manybody/ThreeBodySum.hpp"
#include "jgap/experimental/transform/nbody/2b/CoordinationTransformation.hpp"
#include "jgap/experimental/transform/nbody/3b/MeamTransformation.hpp"

namespace jgap::utils {

    template<size_t Dim, CKernelOfDim<Dim> TKernel>
    std::vector<TwoBodyGapComponent<Dim, TKernel>> createTwoBodyComponents(
        const std::vector<Atoms>& training_data,
        const ValuePtr<TwoBodyTransformation<Dim>>& transformation,
        const TKernel& kernel,
        const Sparsifier<Dim>& sparsifier
    ) {
        std::set<Species2Sorted> all_species_sets;
        Real cutoff = transformation->getCutoffs().maxOverall();

        for (const auto& atoms: training_data) {
            NeighbourLists nl(atoms, cutoff);
            auto sets = Species2Sorted::getAll(nl);
            all_species_sets.insert(sets.begin(), sets.end());
        }

        std::vector<Species2Sorted> species_vec(all_species_sets.begin(), all_species_sets.end());
        std::vector<std::optional<TwoBodyGapComponent<Dim, TKernel>>> components_opt(species_vec.size());

        unseqForIndex(0, species_vec.size(), [&](size_t i) {
            components_opt[i].emplace(species_vec[i], transformation, kernel, sparsifier, training_data);
        });

        std::vector<TwoBodyGapComponent<Dim, TKernel>> components;
        components.reserve(species_vec.size());
        for (auto& comp: components_opt) {
            components.push_back(std::move(*comp));
        }

        return components;
    }

    template<size_t Dim, CKernelOfDim<Dim> TKernel>
    std::vector<ThreeBodyGapComponent<Dim, TKernel>> createThreeBodyComponents(
        const std::vector<Atoms>& training_data,
        const ValuePtr<ThreeBodyTransformation<Dim>>& transformation,
        const TKernel& kernel,
        const Sparsifier<Dim>& sparsifier
    ) {
        std::set<Species3AtomicSorted> all_species_sets;
        Real cutoff = transformation->getCutoffs().maxOverall();

        for (const auto& atoms: training_data) {
            NeighbourLists nl(atoms, cutoff);
            auto sets = Species3AtomicSorted::getAll(nl);
            all_species_sets.insert(sets.begin(), sets.end());
        }

        std::vector<Species3AtomicSorted> species_vec(all_species_sets.begin(), all_species_sets.end());
        std::vector<std::optional<ThreeBodyGapComponent<Dim, TKernel>>> components_opt(species_vec.size());

        unseqForIndex(0, species_vec.size(), [&](size_t i) {
            components_opt[i].emplace(species_vec[i], transformation, kernel, sparsifier, training_data);
        });

        std::vector<ThreeBodyGapComponent<Dim, TKernel>> components;
        components.reserve(species_vec.size());
        for (auto& comp: components_opt) {
            components.push_back(std::move(*comp));
        }

        return components;
    }

    std::map<Species, ValuePtr<NBodyAggregator<1>>> createEamAggregators(
        const ValuePtr<TwoBodyTransformation<1>>& base_pf, const std::vector<Atoms>& training_data, EamMode mode
    );

    template<CKernelOfDim<1> TKernel>
    std::vector<ManyBodyGapComponent<1, TKernel>> createEamComponents(
        const ValuePtr<TwoBodyTransformation<1>>& base_pf,
        const TKernel& kernel,
        const Sparsifier<1>& sparsifier,
        const std::vector<Atoms>& training_data,
        EamMode mode,
        const std::vector<Real>& optional_coeffs = {}
    ) {
        auto aggregators_map = createEamAggregators(base_pf, training_data, mode);
        std::vector<ValuePtr<NBodyAggregator<1>>> aggregators_vec;
        aggregators_vec.reserve(aggregators_map.size());
        for (auto& [species, agg]: aggregators_map) {
            aggregators_vec.push_back(std::move(agg));
        }

        std::vector<std::optional<ManyBodyGapComponent<1, TKernel>>> components_opt(aggregators_vec.size());
        unseqForIndex(0, aggregators_vec.size(), [&](size_t i) {
            components_opt[i].emplace(
                std::move(aggregators_vec[i]), kernel, sparsifier, training_data, optional_coeffs
            );
        });

        std::vector<ManyBodyGapComponent<1, TKernel>> components;
        components.reserve(aggregators_vec.size());
        for (auto& comp: components_opt) {
            components.push_back(std::move(*comp));
        }

        return components;
    }

    template<size_t Dim>
    std::map<Species, ValuePtr<NBodyAggregator<Dim>>> createCoordinationAggregators(
        const ValuePtr<CoordinationTransformation<Dim>>& base_transform, const std::vector<Atoms>& training_data
    ) {
        std::map<Species, ValuePtr<NBodyAggregator<Dim>>> aggregators;

        std::set<Species> all_species;
        for (const auto& atoms: training_data) {
            for (const auto& species: atoms.getSpecies()) {
                all_species.insert(species);
            }
        }

        for (const auto& central_species: all_species) {
            auto aggregator = TwoBodySum<Dim>(central_species);

            for (const auto& contributor_species: all_species) {
                auto transform_clone = base_transform;
                aggregator.extend({central_species, contributor_species}, std::move(transform_clone));
            }
            aggregators.emplace(central_species, std::move(aggregator));
        }

        return aggregators;
    }

    template<size_t Dim, CKernelOfDim<Dim> TKernel>
    std::vector<ManyBodyGapComponent<Dim, TKernel>> createCoordinationComponents(
        const ValuePtr<CoordinationTransformation<Dim>>& base_transform,
        const TKernel& kernel,
        const Sparsifier<Dim>& sparsifier,
        const std::vector<Atoms>& training_data,
        const std::vector<Real>& optional_coeffs = {}
    ) {
        auto aggregators_map = createCoordinationAggregators<Dim>(base_transform, training_data);
        std::vector<ValuePtr<NBodyAggregator<Dim>>> aggregators_vec;
        aggregators_vec.reserve(aggregators_map.size());
        for (auto& [species, agg]: aggregators_map) {
            aggregators_vec.push_back(std::move(agg));
        }

        std::vector<std::optional<ManyBodyGapComponent<Dim, TKernel>>> components_opt(aggregators_vec.size());
        unseqForIndex(0, aggregators_vec.size(), [&](size_t i) {
            components_opt[i].emplace(
                std::move(aggregators_vec[i]), kernel, sparsifier, training_data, optional_coeffs
            );
        });

        std::vector<ManyBodyGapComponent<Dim, TKernel>> components;
        components.reserve(aggregators_vec.size());
        for (auto& comp: components_opt) {
            components.push_back(std::move(*comp));
        }

        return components;
    }

    std::map<Species, ValuePtr<NBodyAggregator<3>>> createMeamAggregators(
        const ValuePtr<MeamTransformation>& base_transform, const std::vector<Atoms>& training_data
    ) {
        std::map<Species, ValuePtr<NBodyAggregator<3>>> aggregators;

        std::set<Species> all_species;
        for (const auto& atoms: training_data) {
            for (const auto& species: atoms.getSpecies()) {
                all_species.insert(species);
            }
        }

        for (const auto& central_species: all_species) {
            auto aggregator = ThreeBodySum<3>(central_species);

            for (const auto& contributor_species1: all_species) {
                for (const auto& contributor_species2: all_species) {
                    auto transform_clone = base_transform;
                    aggregator.extend(
                        {central_species, contributor_species1, contributor_species2}, std::move(transform_clone)
                    );
                }
            }
            aggregators.emplace(central_species, std::move(aggregator));
        }

        return aggregators;
    }

    template<CKernelOfDim<3> TKernel>
    std::vector<ManyBodyGapComponent<3, TKernel>> createMeamComponents(
        const ValuePtr<MeamTransformation>& base_transform,
        const TKernel& kernel,
        const Sparsifier<3>& sparsifier,
        const std::vector<Atoms>& training_data,
        const std::vector<Real>& optional_coeffs = {}
    ) {
        auto aggregators_map = createMeamAggregators(base_transform, training_data);
        std::vector<ValuePtr<NBodyAggregator<3>>> aggregators_vec;
        aggregators_vec.reserve(aggregators_map.size());
        for (auto& [species, agg]: aggregators_map) {
            aggregators_vec.push_back(std::move(agg));
        }

        std::vector<std::optional<ManyBodyGapComponent<3, TKernel>>> components_opt(aggregators_vec.size());
        unseqForIndex(0, aggregators_vec.size(), [&](size_t i) {
            components_opt[i].emplace(
                std::move(aggregators_vec[i]), kernel, sparsifier, training_data, optional_coeffs
            );
        });

        std::vector<ManyBodyGapComponent<3, TKernel>> components;
        components.reserve(aggregators_vec.size());
        for (auto& comp: components_opt) {
            components.push_back(std::move(*comp));
        }

        return components;
    }

}

#endif
