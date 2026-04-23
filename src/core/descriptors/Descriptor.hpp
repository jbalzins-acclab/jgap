#ifndef JGAP_DESCRIPTOR_HPP
#define JGAP_DESCRIPTOR_HPP

#include "../atomic/geometry/Vector3.hpp"
#include "../DataNode.hpp"
#include "core/atomic/energy/Virials.hpp"

#include <span>

namespace jgap {

    template<size_t N_DIMENSIONS>
    struct Descriptor {

        // q
        std::array<double, N_DIMENSIONS> value;

        // dq_i/dH_ab, where H is a matrix s.t.
        // strained cell{\vec{a}, \vec{b}, \vec{c}} = matr{H} * cell{\vec{a}, \vec{b}, \vec{c}}
        std::array<Virials, N_DIMENSIONS> virials;
    };

    template<size_t N_DIMENSIONS, size_t N_ATOMS>
    struct FixedDependenciesDescriptor : Descriptor<N_DIMENSIONS> {

        struct GradientData {
            size_t wrt_atom_index{};
            Vector3 value{};
        };

        std::array<std::array<GradientData, N_ATOMS>, N_DIMENSIONS> gradients{};
    };

    template<size_t N_DIMENSIONS>
    struct ManyBodyDescriptor : Descriptor<N_DIMENSIONS> {
        std::array<std::vector<Vector3>, N_DIMENSIONS> gradients;

        explicit ManyBodyDescriptor(size_t box_n_atoms) {
            for (size_t dim = 0; dim < N_DIMENSIONS; dim++) {
                gradients[dim].resize(box_n_atoms);
            }
        }
    };
}

#endif
