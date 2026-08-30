#ifndef JGAP_STANDARDTABULATION_HPP
#define JGAP_STANDARDTABULATION_HPP

#include <array>
#include <string>
#include "jgap/core/Real.hpp"
#include "jgap/core/potentials/Potential.hpp"
#include "jgap/core/potentials/tabgap/TabGapPotential.hpp"
#include "jgap/core/tabulation/TabulationParams.hpp"
#include "jgap/io/log/CurrentLogger.hpp"
#include "jgap/io/tabgap/TabGapIO.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"
#include "jgap/utils/Utils.hpp"

namespace jgap::utils {

    struct StandardTabulationParams {
        Real r_min_3b = 0.1;
        Real max_eam_density = 10.0;
        size_t n_grid_2b = 5000;
        std::array<size_t, 3> n_grid_3b = {80, 80, 80};
    };

    inline void standardTabulation(
        const Potential& potential, const std::string& output_prefix, const StandardTabulationParams& params = {}
    ) {
        TabulationParams tab_params;
        tab_params.max_cutoffs = potential.getCutoffs();
        tab_params.r_min_3b = params.r_min_3b;
        tab_params.max_eam_density = params.max_eam_density;
        tab_params.n_grid_2b = params.n_grid_2b;
        tab_params.n_grid_3b = params.n_grid_3b;

        const TabGapPotential tabgap{potential.tabulate(tab_params)};
        const Filenames tabgap_files = TabGapIO::write(tabgap, output_prefix);
        JGAP_LOG_INFO("Saved tabulated potential to: {}", vectorToString(tabgap_files));
    }

    inline void standardTabulation(
        const std::string& pot_filename, const std::string& output_prefix, const StandardTabulationParams& params = {}
    ) {
        const ValuePtr<Potential> potential = SerializationRegistry<Potential>::deserialize(pot_filename);
        standardTabulation(*potential, output_prefix, params);
    }
}

#endif
