//
// Created by Jegors Balzins on 18.6.2025.
//

#ifndef TWOBODYKERNELINDEX_HPP
#define TWOBODYKERNELINDEX_HPP

#include <vector>

using namespace std;

namespace jgap {
    struct TwoBodyKernelIndexEntity {
        size_t atomIndex0;
        size_t atomIndex1;

        Vector3 r01;
        double r;

        double fCut;
        double dCut_dr;
    };

    using TwoBodyKernelIndex = vector<TwoBodyKernelIndexEntity>;

    struct TwoBodyDescriptorData {
        double r;
        double fCut;
    };
}

#endif //TWOBODYKERNELINDEX_HPP
