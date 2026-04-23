#ifndef JGAP_TRANSFORMER_HPP
#define JGAP_TRANSFORMER_HPP

#include <vector>

#include "Descriptor.hpp"
#include "core/atomic/geometry/Separations.hpp"
#include "core/atomic/neighbours/NeighbourList.hpp"

namespace jgap {

    class SeparationTransformer {
    public:

        template<size_t N_ATOMS_CONNECTED, size_t TO_DIMENSIONS>
        class FixedDependencies {
        public:
            using TransformFrom = Separations<N_ATOMS_CONNECTED>;
            using TransformTo = FixedDependenciesDescriptor<TO_DIMENSIONS, N_ATOMS_CONNECTED>;

            struct ToBeCalculatedPerDimension {
                Real value;
                std::array<Real, TransformFrom::N_SEPARATIONS> derivatives;
            };

            virtual ~FixedDependencies() = default;

            virtual std::array<ToBeCalculatedPerDimension, TO_DIMENSIONS> evaluate(const TransformFrom& separations)
                = 0;

            TransformTo transform(const TransformFrom& separations);
            std::vector<TransformTo> transform(const std::vector<TransformFrom>& separations);
        };

        template<size_t TO_DIMENSIONS>
        class ManyBody {
        public:
            using TransformTo = std::vector<ManyBodyDescriptor<TO_DIMENSIONS>>;

            virtual ~ManyBody() = default;

            virtual TransformTo transform(const NeighbourList& neighbour_list) = 0;
        };
    };

    template<size_t N_ATOMS_CONNECTED, size_t TO_DIMENSIONS>
    typename SeparationTransformer::FixedDependencies<N_ATOMS_CONNECTED, TO_DIMENSIONS>::TransformTo
    SeparationTransformer::FixedDependencies<N_ATOMS_CONNECTED, TO_DIMENSIONS>::transform(
        const TransformFrom &separations) {

        auto q = evaluate(separations); // q(r_01, r_02, ...) and [dq/dr_ij] : ij in order as in @Separations

        TransformTo result;
        for (size_t dim = 0; dim < TO_DIMENSIONS; dim++) {
            result.value[dim] = q[dim].value;

            auto& grads = result.gradients[dim];
            auto& virials = result.virials[dim];

            for (size_t i = 0; i < TransformFrom::N_ATOMS_CONNECTED; i++) {
                for (size_t j = i + 1; j < TransformFrom::N_ATOMS_CONNECTED; j++) {

                    Vector3& normed_r_ij = separations.between(i, j).direction;

                    grads[j] += normed_r_ij * q[dim].derivatives[flattened_index(i, j)];
                    grads[i] -= normed_r_ij * q[dim].derivatives[flattened_index(i, j)];

                    virials += q[dim].derivatives[flattened_index(i, j)]
                            * separations.between(i, j).virials;
                }
            }
        }

        return result;
    }

    template<size_t N_ATOMS_CONNECTED, size_t TO_DIMENSIONS>
    std::vector<typename SeparationTransformer::FixedDependencies<N_ATOMS_CONNECTED, TO_DIMENSIONS>::TransformTo>
    SeparationTransformer::FixedDependencies<N_ATOMS_CONNECTED, TO_DIMENSIONS>::transform(
        const std::vector<TransformFrom> &separations) {
        return std::ranges::transform(
            separations.begin(), separations.end(),
            [this](const TransformFrom& from) -> auto { return transform(from); }
        );
    }
}


#endif //JGAP_TRANSFORMER_HPP
