#ifndef JGAP_TABULATIONPARAMS_HPP
#define JGAP_TABULATIONPARAMS_HPP
#include <set>

#include "core/atomic/species/Species.hpp"
#include "core/potentials/Cutoffs.hpp"

namespace jgap {
    struct TabulationParams {
        Cutoffs max_cutoffs;
        Real max_eam_density = 10.0;
        std::set<Species> species;
    };
}

#endif
