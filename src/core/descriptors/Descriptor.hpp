#ifndef DESCRIPTOR_HPP
#define DESCRIPTOR_HPP

#include "data/Vector3.hpp"
#include "memory/MatrixBlock.hpp"

#include <memory>
#include <vector>
#include <string>
#include "data/DataNode.hpp"

#include "../kernels/Kernel.hpp"
#include "data/CutoffRanges.hpp"
#include "../../data/tabulation/TabulationData.hpp"

namespace jgap {

    template<size_t N_DIMENSIONS, size_t N_GRADIENTS>
    struct Descriptor {
        struct GradientData {
            size_t wrt_atom_index{};
            Vector3 value{};
        };

        std::array<double, N_DIMENSIONS> value; // q
        std::array<std::array<GradientData, N_DIMENSIONS>, N_GRADIENTS> gradients; // dq_i/dr_i
        // dq_i/dH_ab, where H is a matrix s.t.
        // strained cell{\vec{a}, \vec{b}, \vec{c}} = matr{H} * cell{\vec{a}, \vec{b}, \vec{c}}
        std::array<Virials, N_GRADIENTS> virials;
    };
    /*template<size_t N_DIMENSIONS, size_t N_GRADIENTS>
        class DescriptorFinder {
        public:
            using DescriptorsFiltered = std::map<
                EncodedSpeciesSets,
                Descriptor<N_DIMENSIONS, N_GRADIENTS>
            >;

            virtual ~DescriptorFinder() = default;
            virtual DescriptorsFiltered find(const AtomicStructure &atomic_structure) = 0;
        };*/
    /*
    template<size_t N_DIMENSIONS, size_t N_STATIC_GRADIENT_SLOTS>
    struct DynamicDescriptor {
        struct GradientData {
            Vector3 value;
            size_t wrt_atom_index;
        };

        std::array<double, N_DIMENSIONS> value; // q

        size_t n_gradients;
        std::array<std::array<GradientData, N_DIMENSIONS>, N_STATIC_GRADIENT_SLOTS> gradients_static;
        std::vector<std::array<GradientData, N_DIMENSIONS>> gradients_dynamic;
    };*/
}

#endif
