#ifndef JGAP_NEIGHBOURDATA_HPP
#define JGAP_NEIGHBOURDATA_HPP

#include <array>
#include <vector>
#include <optional>
#include <numeric>
#include <stdexcept>
#include "../geometry/Separations.hpp"
#include "../geometry/Vector3.hpp"

namespace jgap {
    struct NeighbourData {
        size_t atom_index;
        Separation separation;
    };
}

#endif