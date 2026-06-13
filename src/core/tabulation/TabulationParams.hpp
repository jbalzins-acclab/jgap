#ifndef JGAP_TABULATIONPARAMS_HPP
#define JGAP_TABULATIONPARAMS_HPP

#include <set>

#include "core/atomic/species/Species.hpp"
#include "core/potentials/Cutoffs.hpp"

namespace jgap {
    struct TabulationParams {
        Cutoffs max_cutoffs;
        Real max_eam_density{12.0};

        size_t n_grid_2b{5000};
        std::array<size_t, 3> n_grid_3b{100, 100, 100};
    };
}

#endif
