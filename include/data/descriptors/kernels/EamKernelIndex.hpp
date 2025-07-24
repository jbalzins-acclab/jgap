#ifndef EAMKERNELINDEX_HPP
#define EAMKERNELINDEX_HPP

#include <vector>
#include <map>

#include "data/BasicDataTypes.hpp"

using namespace std;

namespace jgap {
    struct EamDensityData {
        // to avoid overhead vector.index =
        size_t atAtomIndex;
        double density; // rho_i = sum(...)
        vector<pair<NeighbourData, double>> densityDerivatives; // drho_i / dr_ij (j = NeighbourData.index)
    };

    using EamKernelIndexPerSpecies = vector<EamDensityData>;
    using EamKernelIndex = map<Species, EamKernelIndexPerSpecies>;
}

#endif //EAMKERNELINDEX_HPP
