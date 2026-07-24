#ifndef JGAP_CLUSTER2_HPP
#define JGAP_CLUSTER2_HPP

#include <array>
#include <cstddef>
#include "jgap/core/Real.hpp"

namespace jgap {
    struct Cluster2 {
        std::array<size_t, 2> atom_indexes;
        Real r01;
    };
}

#endif
