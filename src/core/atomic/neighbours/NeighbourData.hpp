#ifndef JGAP_NEIGHBOURDATA_HPP
#define JGAP_NEIGHBOURDATA_HPP

#include "../geometry/Separation.hpp"

namespace jgap {
    struct NeighbourData {
        size_t atom_index{};
        Separation separation{};
    };
}

#endif