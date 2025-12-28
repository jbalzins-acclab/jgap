#ifndef EAMKERNELINDEX_HPP
#define EAMKERNELINDEX_HPP

#include <vector>
#include <map>

#include "data/Vector3.hpp"

namespace jgap {
    struct EamDensityData {
        // to avoid overhead std::vector.index =
        size_t atAtomIndex;
        double density; // rho_i = sum(...)
        std::vector<std::pair<NeighbourData, double>> densityDerivatives; // drho_i / dr_ij (j = NeighbourData.index)
    };

    using EamKernelIndexPerSpecies = std::vector<EamDensityData>;
    using EamKernelIndex = std::map<Species, EamKernelIndexPerSpecies>;
}

#endif
