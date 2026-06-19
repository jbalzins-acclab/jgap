#ifndef JGAP_STANDARDGAPPARAMS_HPP
#define JGAP_STANDARDGAPPARAMS_HPP

#include <array>
#include <optional>
#include <string>
#include "core/transform/eam/EamPairFunction.hpp"
#include "core/transform/eam/FSGenPairFunction.hpp"
#include "core/fit/gap/regularization/SimpleRegularizationRules.hpp"

namespace jgap {
    struct StandardGapParams {
        size_t seed;
        // External ZBL: read its coefficients from this file if set, otherwise use the built-in dataset.
        std::optional<std::string> zbl_dataset_file{};
        // EAM params
        ValuePtr<ClusterTransformation<1, 2>> eam_pf = FSGenPairFunction(4.5, 3.0);
        size_t eam_n_sparse = 20;
        EamMode eam_mode = EamMode::Blind;
        // 3-body params
        Real cutoff3 = 3.7;
        Real cutoff3_width = 0.6;
        size_t n_sparse3 = 500;
        std::array<bool, 4> active3 = {true, true, true, false};
        // 2-body params
        Real cutoff2 = 4.5;
        Real cutoff2_width = 1.0;
        size_t n_sparse2 = 20;
        std::array<bool, 2> active2 = {true, false};
        // Regularization
        ValuePtr<RegularizationRules> regularization_rules = SimpleRegularizationRules();
    };
}

#endif