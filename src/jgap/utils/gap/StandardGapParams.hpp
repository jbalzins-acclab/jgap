#ifndef JGAP_STANDARDGAPPARAMS_HPP
#define JGAP_STANDARDGAPPARAMS_HPP

#include <array>
#include <optional>
#include <string>
#include <vector>
#include "jgap/core/Real.hpp"
#include "jgap/core/atomic/species/Species.hpp"
#include "jgap/core/transform/nbody/2b/eam/EamPairFunction.hpp"

namespace jgap::utils {

    enum class EamPairFunctionType { FSGen2, FSGen3, Coscutoff, Polycutoff };

    struct StandardGapParams {
        size_t seed = 42;

        // External ScreenedCoulomb: read its coefficients from this file if set, otherwise use the built-in dataset.
        std::optional<std::string> screened_coulomb_dataset_file{};

        // 2-body params
        Real cutoff2 = 4.5;
        Real cutoff2_width = 1.0;
        size_t n_sparse2 = 20;

        // EAM params
        EamMode eam_mode = EamMode::Blind;
        EamPairFunctionType eam_pair_function = EamPairFunctionType::FSGen3;
        size_t eam_n_sparse = 20;
        Real eam_min_density = 0.05;

        // 3-body params
        Real cutoff3 = 3.7;
        Real cutoff3_width = 0.6;
        size_t n_sparse3 = 500;

        // Mandatory RAM limit in gigabytes for ElementalQRGapFit out-of-core execution.
        Real approx_ram_limit_gb = 4.0;
    };
}

#endif
