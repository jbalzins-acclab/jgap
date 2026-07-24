#ifndef JGAP_CLUSTER3_HPP
#define JGAP_CLUSTER3_HPP

#include <array>
#include <cassert>
#include <cstddef>
#include "jgap/core/Real.hpp"

#include "Cluster.hpp"

namespace jgap {

    struct Cluster3 {
        std::array<size_t, 3> atom_indexes;
        std::array<Real, 3> separation_magnitudes;

        Cluster3() = default;
        Cluster3(const Cluster3 &) = default;

        Real separationBetween(const size_t lower_index, const size_t higher_index) const {
            assert(lower_index < higher_index);
            return separation_magnitudes[flattenedIndex(lower_index, higher_index)];
        }

        Real &separationBetween(const size_t lower_index, const size_t higher_index) {
            assert(lower_index < higher_index);
            return separation_magnitudes[flattenedIndex(lower_index, higher_index)];
        }
    };
}

#endif
