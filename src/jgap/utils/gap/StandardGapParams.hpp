#ifndef JGAP_STANDARDGAPPARAMS_HPP
#define JGAP_STANDARDGAPPARAMS_HPP

#include <array>
#include <optional>
#include <string>
#include <vector>
#include "jgap/core/atomic/species/Species.hpp"
#include "jgap/core/atomic/species/composition/Species2Sorted.hpp"
#include "jgap/core/atomic/species/composition/Species3AtomicSorted.hpp"
#include "jgap/core/transform/nbody/2b/eam/EamPairFunction.hpp"

namespace jgap::utils {

    enum class EamPairFunctionType { FSGen2, FSGen3, Coscutoff, Polycutoff };

    enum class ThreeBodyTransformationType { Angle, Distances };

    struct StandardGap2bParams {
        std::optional<Species2Sorted> species{}; // If not provided, applies as non-species default
        double cutoff = 4.5;
        double cutoff_width = 1.0;
        size_t n_sparse = 20;
        double energy_scale = 10.0;
        double length_scale = 1.0;

        bool operator==(const StandardGap2bParams& other) const = default;
    };

    struct StandardGapEamParams {
        std::optional<Species> species{}; // If not provided, applies as non-species default
        EamMode eam_mode = EamMode::Blind;
        EamPairFunctionType eam_pair_function = EamPairFunctionType::FSGen3;
        double cutoff = 4.5;
        size_t n_sparse = 20;
        double min_density = 0.05;
        double energy_scale = 1.0;
        double length_scale = 1.0;

        bool operator==(const StandardGapEamParams& other) const = default;
    };

    struct StandardGap3bParams {
        std::optional<Species3AtomicSorted> species{}; // If not provided, applies as non-species default
        ThreeBodyTransformationType transformation_type = ThreeBodyTransformationType::Angle;
        double cutoff = 3.7;
        double cutoff_width = 0.6;
        size_t n_sparse = 500;
        double energy_scale = 1.0;
        std::array<double, 3> length_scales = {1.0, 1.0, 1.0};

        bool operator==(const StandardGap3bParams& other) const = default;
    };

    struct StandardGapParams {
        size_t seed = 42;

        // External ScreenedCoulomb: read its coefficients from this file if set, otherwise use the built-in dataset.
        std::optional<std::string> screened_coulomb_dataset_file{};

        // Mandatory RAM limit in gigabytes for ElementIncrementalQRGapFit out-of-core execution.
        double approx_ram_limit_gb = 4.0;

        // Default non-species parameters (applied to all combos without specific params).
        // Can be set to std::nullopt to erase / disable default components.
        std::optional<StandardGap2bParams> default_2b = StandardGap2bParams{};
        std::optional<StandardGapEamParams> default_eam = StandardGapEamParams{};
        std::optional<StandardGap3bParams> default_3b = StandardGap3bParams{};

        // Species-specific parameters (multiple per species combo allowed)
        std::vector<StandardGap2bParams> species_2b = {};
        std::vector<StandardGapEamParams> species_eam = {};
        std::vector<StandardGap3bParams> species_3b = {};
    };
}

#endif
