#ifndef JGAP_SEPARATIONS_HPP
#define JGAP_SEPARATIONS_HPP

#include "Vector3.hpp"
#include "Separation.hpp"
#include <cassert>

namespace jgap {

    static constexpr size_t flattened_index(size_t lower_index, size_t higher_index) {
#ifdef DEBUG
        assert(lower_index < higher_index);
#endif
        // TODO: debug
        return higher_index*(higher_index-1)/2 + lower_index;
    }

    template<size_t N_ATOMS_CONNECTED>
    requires(N_ATOMS_CONNECTED > 1)
    struct Separations {
        static constexpr size_t N_SEPARATIONS = N_ATOMS_CONNECTED * (N_ATOMS_CONNECTED - 1) / 2;

        std::array<size_t, N_ATOMS_CONNECTED> atom_indexes;
        std::array<Separation, N_SEPARATIONS> separations;

        const Separation& between(const size_t lower_index, const size_t higher_index) const {
            return separations[flattened_index(lower_index, higher_index)];
        }

        Separation& between(const size_t lower_index, const size_t higher_index) {
            return separations[flattened_index(lower_index, higher_index)];
        }
    };
}

#endif