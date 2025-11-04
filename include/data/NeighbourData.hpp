#ifndef JGAP_NEIGHBOURDATA_HPP
#define JGAP_NEIGHBOURDATA_HPP
#include "Vector3.hpp"

namespace jgap {
    struct NeighbourData {
        size_t index;
        Vector3 offset;

        double distance;
    };

    using NeighboursData = vector<NeighbourData>;
}

#endif