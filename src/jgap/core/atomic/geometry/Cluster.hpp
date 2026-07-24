#ifndef JGAP_CLUSTER_HPP
#define JGAP_CLUSTER_HPP

#include <cassert>
#include <optional>
#include <vector>

namespace jgap {

    /// Calculates flat-array index for a symmetric matrix.
    ///
    /// \note The first matrix index must be lower than the second index.
    ///      To save performance, this is unchecked unless the DEBUG flag is on.
    static constexpr size_t flattenedIndex(size_t lower_index, size_t higher_index) {
        assert(lower_index < higher_index);
        return higher_index * (higher_index - 1) / 2 + lower_index;
    }

    static constexpr size_t nSeparations(size_t cluster_size) { return cluster_size * (cluster_size - 1) / 2; }
}

#endif
