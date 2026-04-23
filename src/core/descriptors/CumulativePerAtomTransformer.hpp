#ifndef JGAP_CUMULATIVETRANSFORMER_HPP
#define JGAP_CUMULATIVETRANSFORMER_HPP
#include "SeparationTransformer.hpp"

namespace jgap {

    template<size_t N_ATOMS_CONNECTED, size_t TO_DIMENSIONS>
    class CumulativePerAtomTransformer : public SeparationTransformer::ManyBody<TO_DIMENSIONS> {
    public:
        using TransformTo = SeparationTransformer::ManyBody<TO_DIMENSIONS>::TransformTo;

        struct ToBeCalculatedPerDimension {
            Real value;
            std::array<Real, Separations<N_ATOMS_CONNECTED>::N_SEPARATIONS> derivatives;
        };

        ~CumulativePerAtomTransformer() override = default;

        TransformTo transform(const NeighbourList& neighbour_list) override {
            TransformTo result(neighbour_list.nAtoms());

            for (const Separations<N_ATOMS_CONNECTED>& sep: separations) {
                // q(r_01, r_02, ...) and [dq/dr_ij] : ij in order as in @Separations
                auto q = contribution(sep);

                for (size_t dim = 0; dim < TO_DIMENSIONS; dim++) {
                    result.value[dim] += q[dim].value;


                }
            }

            return result;
        }

        virtual std::array<ToBeCalculatedPerDimension, TO_DIMENSIONS> contribution(
            const Separations<N_ATOMS_CONNECTED>& sep);
    };
}

#endif
