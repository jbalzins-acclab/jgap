#ifndef JGAP_STANDARDGAPFIT_HPP
#define JGAP_STANDARDGAPFIT_HPP

#include <set>
#include <utility>
#include "StandardGapParams.hpp"
#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/impl/cutoff/CosCutoff.hpp"
#include "jgap/impl/kernels/SquaredExpKernel.hpp"
#include "jgap/core/potentials/CompositePotential.hpp"
#include "jgap/core/potentials/coulomb/ScreenedCoulombPotential.hpp"
#include "jgap/core/potentials/gap/GapPotential.hpp"
#include "jgap/core/potentials/gap/component/ManyBodyGapComponent.hpp"
#include "jgap/core/potentials/gap/component/ThreeBodyGapComponent.hpp"
#include "jgap/core/potentials/gap/component/TwoBodyGapComponent.hpp"
#include "jgap/core/potentials/isolated/IsolatedAtomPotential.hpp"
#include "jgap/core/sparsification/HistogramUniformSparsifier.hpp"
#include "jgap/impl/transform/nbody/2b/PairDistanceTransformation.hpp"
#include "jgap/impl/transform/nbody/2b/eam/CoscutoffPairFunction.hpp"
#include "jgap/impl/transform/nbody/2b/eam/FSGenPairFunction.hpp"
#include "jgap/impl/transform/nbody/2b/eam/PolycutoffPairFunction.hpp"
#include "jgap/impl/transform/nbody/3b/Angle3bTransformation.hpp"
#include "jgap/impl/transform/nbody/3b/Distances3bTransformation.hpp"
#include "jgap/impl/transform/manybody/TwoBodySum.hpp"
#include "jgap/impl/fit/gap/ElementIncrementalQRGapFit.hpp"
#include "jgap/impl/fit/gap/BlockIncrementalQRGapFit.hpp"
#include "jgap/impl/fit/gap/QRGapFit.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"
#include "jgap/utils/gap/GapComponentUtils.hpp"

namespace jgap::utils {

    inline ValuePtr<EamPairFunction> makeStandardEamPairFunction(EamPairFunctionType type, double cutoff) {
        switch (type) {
            case EamPairFunctionType::FSGen2:
                return FSGenPairFunction(cutoff, 2.0);
            case EamPairFunctionType::FSGen3:
                return FSGenPairFunction(cutoff, 3.0);
            case EamPairFunctionType::Coscutoff:
                return CoscutoffPairFunction(cutoff, 0.0);
            case EamPairFunctionType::Polycutoff:
                return PolycutoffPairFunction(cutoff, 0.0);
        }
        std::unreachable();
    }

    inline ValuePtr<ThreeBodyTransformation<4>> makeStandard3bTransformation(
        ThreeBodyTransformationType type,
        double cutoff,
        double cutoff_width
    ) {
        auto cut = CosCutoff(cutoff, cutoff_width);
        switch (type) {
            case ThreeBodyTransformationType::Angle:
                return Angle3bTransformation(cut);
            case ThreeBodyTransformationType::Distances:
                return Distances3bTransformation(cut);
        }
        std::unreachable();
    }

    inline void standardGapFit(
        const std::string& filename,
        const std::vector<Atoms>& training_data,
        const std::vector<Regularization>& sigmas,
        const StandardGapParams& params = {}
    ) {
        if (training_data.empty()) {
            JGAP_LOG_AND_THROW("Training data cannot be empty");
        }

        GapPotential potential;

        std::set<Species> all_species;
        for (const auto& atoms: training_data) {
            for (const auto& s: atoms.getSpecies()) {
                all_species.insert(s);
            }
        }

        // ====================================================================================
        // 2-Body Components
        // ====================================================================================
        std::set<Species2Sorted> specific_2b_pairs;
        for (const auto& p: params.species_2b) {
            if (p.species.has_value()) {
                specific_2b_pairs.insert(*p.species);
            }
        }

        // Add species-specific 2-body components
        for (const auto& p: params.species_2b) {
            if (!p.species.has_value() || p.n_sparse == 0) continue;
            const auto& pair = *p.species;

            auto trans2 = PairDistanceTransformation(CosCutoff(p.cutoff, p.cutoff_width));
            auto kernel2 = SquaredExpKernel<1, 1>(p.energy_scale, {p.length_scale});
            auto sparsifier2 = HistogramUniformSparsifier<2>(params.seed, p.n_sparse, std::array{true, false});
            TwoBodyGapComponent<2, SquaredExpKernel<1, 1>> comp(
                pair, trans2, kernel2, sparsifier2, training_data
            );
            if (comp.nSparsePoints() == 0) continue;

            potential.addComponent(std::move(comp));
        }

        // Add default 2-body components for all pairs NOT in specific_2b_pairs
        if (params.default_2b.has_value() && params.default_2b->n_sparse > 0) {
            const auto& def_p = *params.default_2b;
            std::set<Species2Sorted> candidate_pairs;
            for (const auto& atoms: training_data) {
                NeighbourLists nl(atoms, def_p.cutoff);
                auto sets = Species2Sorted::getAll(nl);
                candidate_pairs.insert(sets.begin(), sets.end());
            }

            for (const auto& pair: candidate_pairs) {
                if (specific_2b_pairs.contains(pair)) {
                    continue; // Already handled by specific params
                }

                auto trans2 = PairDistanceTransformation(CosCutoff(def_p.cutoff, def_p.cutoff_width));
                auto kernel2 = SquaredExpKernel<1, 1>(def_p.energy_scale, {def_p.length_scale});
                auto sparsifier2 = HistogramUniformSparsifier<2>(params.seed, def_p.n_sparse, std::array{true, false});
                TwoBodyGapComponent<2, SquaredExpKernel<1, 1>> comp(
                    pair, trans2, kernel2, sparsifier2, training_data
                );
                if (comp.nSparsePoints() == 0) continue;

                potential.addComponent(std::move(comp));
            }
        }

        // ====================================================================================
        // ManyBodyGapComponent with EAM Pair Function
        // ====================================================================================
        std::set<Species> specific_eam_species;
        for (const auto& p: params.species_eam) {
            if (p.species.has_value()) {
                specific_eam_species.insert(*p.species);
            }
        }

        auto buildEamComponent = [&](const Species& central_species, const StandardGapEamParams& p) {
            if (p.n_sparse == 0) return;
            auto eam_pf = makeStandardEamPairFunction(p.eam_pair_function, p.cutoff);

            auto aggregator = TwoBodySum<1>(central_species);
            auto Z_center_opt = central_species.atomicNumber();
            if (!Z_center_opt && p.eam_mode != EamMode::Blind && p.eam_mode != EamMode::EAM) {
                JGAP_LOG_AND_THROW("Central species of unknown atomic number - incompatible with the EAM mode");
            }
            double Z_center = static_cast<double>(Z_center_opt.value_or(0));

            for (const auto& contributor_species: all_species) {
                auto pf_clone = eam_pf;
                auto& eam_pf_clone = dynamic_cast<EamPairFunction&>(*pf_clone);

                auto Z_contrib_opt = contributor_species.atomicNumber();
                if (!Z_contrib_opt && p.eam_mode != EamMode::Blind) {
                    JGAP_LOG_AND_THROW("Contributor species of unknown atomic number - incompatible with the EAM mode");
                }

                double prefactor = 1.0;
                double Z_contrib = static_cast<double>(Z_contrib_opt.value_or(0));
                if (p.eam_mode == EamMode::FSsym) {
                    prefactor = std::sqrt(Z_contrib * Z_center) / 40.0;
                } else if (p.eam_mode == EamMode::FSgen) {
                    prefactor = std::pow(Z_center, 0.1) * std::sqrt(Z_contrib) / 10.0;
                } else if (p.eam_mode == EamMode::EAM) {
                    prefactor = std::sqrt(Z_contrib) / 10.0;
                }

                eam_pf_clone.setPrefactor(prefactor);
                aggregator.extend({central_species, contributor_species}, std::move(pf_clone));
            }

            auto kernel_eam = SquaredExpKernel<1, 0>(p.energy_scale, {p.length_scale});
            auto sparsifier_eam = HistogramUniformSparsifier<1>(
                params.seed, p.n_sparse, std::nullopt, std::nullopt, Descriptor<1>{p.min_density}
            );

            ValuePtr<NBodyAggregator<1>> agg_ptr = std::move(aggregator);
            ManyBodyGapComponent<1, SquaredExpKernel<1, 0>> comp(
                agg_ptr, kernel_eam, sparsifier_eam, training_data
            );
            if (comp.nSparsePoints() == 0) return;

            potential.addComponent(std::move(comp));
        };

        // Add species-specific EAM components
        for (const auto& p: params.species_eam) {
            if (p.species.has_value() && all_species.contains(*p.species)) {
                buildEamComponent(*p.species, p);
            }
        }

        // Add default EAM components for species not in specific_eam_species
        if (params.default_eam.has_value()) {
            for (const auto& central_species: all_species) {
                if (!specific_eam_species.contains(central_species)) {
                    buildEamComponent(central_species, *params.default_eam);
                }
            }
        }

        // ====================================================================================
        // 3-Body Components
        // ====================================================================================
        std::set<Species3AtomicSorted> specific_3b_triplets;
        for (const auto& p: params.species_3b) {
            if (p.species.has_value()) {
                specific_3b_triplets.insert(*p.species);
            }
        }

        // Add species-specific 3-body components
        for (const auto& p: params.species_3b) {
            if (!p.species.has_value() || p.n_sparse == 0) continue;
            const auto& triplet = *p.species;

            auto trans3 = makeStandard3bTransformation(p.transformation_type, p.cutoff, p.cutoff_width);
            auto kernel3 = SquaredExpKernel<3, 1>(p.energy_scale, p.length_scales);
            auto sparsifier3 =
                HistogramUniformSparsifier<4>(params.seed, p.n_sparse, std::array{true, true, true, false});
            ThreeBodyGapComponent<4, SquaredExpKernel<3, 1>> comp(
                triplet, trans3, kernel3, sparsifier3, training_data
            );
            if (comp.nSparsePoints() == 0) continue;

            potential.addComponent(std::move(comp));
        }

        // Add default 3-body components for triplets not in specific_3b_triplets
        if (params.default_3b.has_value() && params.default_3b->n_sparse > 0) {
            const auto& def_p = *params.default_3b;
            std::set<Species3AtomicSorted> candidate_triplets;
            for (const auto& atoms: training_data) {
                NeighbourLists nl(atoms, def_p.cutoff);
                auto sets = Species3AtomicSorted::getAll(nl);
                candidate_triplets.insert(sets.begin(), sets.end());
            }

            for (const auto& triplet: candidate_triplets) {
                if (specific_3b_triplets.contains(triplet)) {
                    continue; // Already handled by specific params
                }

                auto trans3 = makeStandard3bTransformation(def_p.transformation_type, def_p.cutoff, def_p.cutoff_width);
                auto kernel3 = SquaredExpKernel<3, 1>(def_p.energy_scale, def_p.length_scales);
                auto sparsifier3 =
                    HistogramUniformSparsifier<4>(params.seed, def_p.n_sparse, std::array{true, true, true, false});
                ThreeBodyGapComponent<4, SquaredExpKernel<3, 1>> comp(
                    triplet, trans3, kernel3, sparsifier3, training_data
                );
                if (comp.nSparsePoints() == 0) continue;

                potential.addComponent(std::move(comp));
            }
        }

        if (potential.getComponents().empty()) {
            JGAP_LOG_AND_THROW("Cannot make a standard GAP potential without any components");
        }

        IsolatedAtomPotential isolated_atom_pot{training_data};
        ScreenedCoulombPotential sc_pot =
            params.screened_coulomb_dataset_file
                ? ScreenedCoulombPotential{*params.screened_coulomb_dataset_file, training_data}
                : ScreenedCoulombPotential{training_data};

        CompositePotential external{{
            {"isolated", isolated_atom_pot},
            {"screened_coulomb", sc_pot},
        }};

        potential.optional_external_potential = external;

        // Always perform ElementIncrementalQRGapFit with single-element structures first
        ElementIncrementalQRGapFit fitter(1e-8, params.approx_ram_limit_gb);
        fitter.fit(potential, training_data, sigmas);

        SerializationRegistry<Potential>::serialize(potential, filename);
        JGAP_LOG_INFO("Saved fitted potential to {}", filename);
    }
}

#endif
