#ifndef JGAP_CLUSTER3_HPP
#define JGAP_CLUSTER3_HPP

#include <array>
#include <cassert>
#include <concepts>
#include <cstddef>
#include "jgap/core/atomic/geometry/Separation.hpp"

namespace jgap {
    struct Cluster3 {
        std::array<size_t, 3> atom_indexes;
        std::array<Separation, 3> separations;

        Separation& separation01() { return separations[0]; }
        Separation& separation02() { return separations[1]; }
        Separation& separation12() { return separations[2]; }

        const Separation& separation01() const { return separations[0]; }
        const Separation& separation02() const { return separations[1]; }
        const Separation& separation12() const { return separations[2]; }

        Separation& separation(size_t flat_index) { return separations[flat_index]; }
        const Separation& separation(size_t flat_index) const { return separations[flat_index]; }

        std::array<size_t, 2> atomIndexes(size_t flat_index) const {
            static constexpr std::array<std::array<size_t, 2>, 3> index_pairs = {{{0, 1}, {0, 2}, {1, 2}}};
            const auto [i, j] = index_pairs[flat_index];
            return {atom_indexes[i], atom_indexes[j]};
        }
    };

    template<typename Func>
    concept Cluster3Callback = std::invocable<Func, const Cluster3&>;
}

#endif
