//
// Created by Jegors Balzins on 19.6.2025.
//

#ifndef EAMKERNELINDEX_HPP
#define EAMKERNELINDEX_HPP

#include <memory>
#include <vector>

using namespace std;

namespace jgap {
    struct EamDensityData {
        // to avoid overhead vector.index =
        double density; // rho_i = sum(...)
        vector<pair<NeighbourData, double>> densityDerivatives; // drho_i / dr_ij (j = NeighbourData.index)
    };

     using EamKernelIndex = vector<EamDensityData>;
}

#endif //EAMKERNELINDEX_HPP
