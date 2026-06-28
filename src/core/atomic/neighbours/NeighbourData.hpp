#ifndef JGAP_NEIGHBOURDATA_HPP
#define JGAP_NEIGHBOURDATA_HPP

#include "../geometry/Separation.hpp"

namespace jgap {
    struct NeighbourData {
        /// Index corresponds to the underlying \Atoms object.
        size_t neighbour_index{};
        Separation separation{};
    };
}

#endif