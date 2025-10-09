#ifndef TWOBODYINDEX_HPP
#define TWOBODYINDEX_HPP

#include <utility>
#include <vector>

#include "data/BasicDataTypes.hpp"

using namespace std;

namespace jgap {
    struct TwoBodyIndexEntity {
        size_t atomIndex0;
        size_t atomIndex1;

        Vector3 r01;
        double r;

        double fCut;
        double dCut_dr;
    };

    using TwoBodyIndex = vector<TwoBodyIndexEntity>;

    struct TwoBodyDescriptorData {
        double r;
        double fCut;
    };
}

#endif
