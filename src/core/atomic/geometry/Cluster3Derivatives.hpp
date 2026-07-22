#ifndef JGAP_CLUSTER3DERIVATIVES_HPP
#define JGAP_CLUSTER3DERIVATIVES_HPP

#include <array>
#include <cassert>
#include <cstddef>
#include "Separation.hpp"

#include "Cluster.hpp"

namespace jgap {

    struct Cluster3Derivatives {
        std::array<SeparationDerivatives, 3> val;

        const SeparationDerivatives &derivativesBetween(const size_t lower_index, const size_t higher_index) const {
            assert(lower_index < higher_index);
            return val[flattenedIndex(lower_index, higher_index)];
        }

        SeparationDerivatives &derivativesBetween(const size_t lower_index, const size_t higher_index) {
            assert(lower_index < higher_index);
            return val[flattenedIndex(lower_index, higher_index)];
        }
    };
}

#endif
