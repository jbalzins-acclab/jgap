#ifndef JGAP_STANDARDGAPPARAMS_HPP
#define JGAP_STANDARDGAPPARAMS_HPP

#include <array>
#include <optional>
#include <string>
#include "jgap/core/ValuePtr.hpp"
#include "jgap/core/fit/gap/regularization/SimpleRegularizationRules.hpp"
#include "jgap/core/transform/nbody/2b/eam/EamPairFunction.hpp"
#include "jgap/core/transform/nbody/2b/eam/FSGenPairFunction.hpp"

namespace jgap::utils {
    struct StandardGapParams {

        size_t seed;

        // External ScreenedCoulomb: read its coefficients from this file if set, otherwise use the built-in dataset.
        std::optional<std::string> screened_coulomb_dataset_file{};

        // 2-body params
        Real cutoff2 = 4.5;
        Real cutoff2_width = 1.0;
        size_t n_sparse2 = 20;

        // EAM params
        EamMode eam_mode = EamMode::Blind;
        ValuePtr<EamPairFunction> eam_pf = FSGenPairFunction(4.5, 3.0);
        size_t eam_n_sparse = 20;
        Real eam_min_density = 0.05;

        // 3-body params
        Real cutoff3 = 3.7;
        Real cutoff3_width = 0.6;
        size_t n_sparse3 = 500;

        // Regularization
        ValuePtr<RegularizationRules> regularization_rules = SimpleRegularizationRules();
    };
}

#endif
