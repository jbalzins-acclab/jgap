#ifndef JGAP_CLUSTER_HPP
#define JGAP_CLUSTER_HPP

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

    template<size_t NAtoms>
    requires(NAtoms > 1)
    struct Cluster {
        static constexpr size_t NSeparations = NAtoms * (NAtoms - 1) / 2;

        std::array<size_t, NAtoms> atom_indexes;
        std::array<Separation, NSeparations> separations;

        const Separation& between(const size_t lower_index, const size_t higher_index) const {
            return separations[flattened_index(lower_index, higher_index)];
        }

        Separation& between(const size_t lower_index, const size_t higher_index) {
            return separations[flattened_index(lower_index, higher_index)];
        }
    };
}

#endif