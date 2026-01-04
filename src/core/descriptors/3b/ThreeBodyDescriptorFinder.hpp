#ifndef THREEBODYDESCRIPTOR_HPP
#define THREEBODYDESCRIPTOR_HPP

#include <queue>
#include <utility>

#include "core/descriptors/Descriptor.hpp"
#include "../../kernels/3b/ThreeBodySE.hpp"
#include "core/cutoff/CosCutoff.hpp"
#include "core/descriptors/DescriptorFinder.hpp"
#include "core/kernels/KernelCollection.hpp"
#include "core/sparsification/Sparsifier.hpp"
#include "core/tabulation/Tabulatable.hpp"
#include "data/descriptors/kernels/ThreeBodyIndex.hpp"

#include "io/parse/ParserRegistry.hpp"
#include "memory/MatrixBlock.hpp"

namespace jgap {

    /**
     * val = (r_01, r_02, r_12),
     * gradients = (
     *  {0, { grad_r0(r_01),  grad_r0(r_02),  grad_r0(r_12) }},
     *  {0, { grad_r0(r_01),  grad_r0(r_02),  grad_r0(r_12) }},
     *  {0, { grad_r0(r_01),  grad_r0(r_02),  grad_r0(r_12) }}
     * ),
     * virials = (
     * )
     */
    using NewThreeBodyDescriptor = Descriptor<3, 3>;
    using ThreeBodyDescriptorsFiltered = std::map<EncodedSpeciesSets, std::vector<NewThreeBodyDescriptor>>;

    struct ThreeBodySeparationsTransformed {
        std::array<double, 3> q{};
        std::array<Vector3, 3> dq_dr_ij{};
    };

    class ThreeBodyDescriptorFinder/*, Serializable, Tabulatable*/ {
    public:
        //SETUP_PARSER_AND_SERIALIZATION(ThreeBodyDescriptorFinder, ThreeBodyDescriptorFinder, 3b)

        std::shared_ptr<CutoffFunction> cutoff_function;
        bool calculate_derivatives;

        ThreeBodyDescriptorFinder(std::shared_ptr<CutoffFunction> cutoff_function, bool calculate_derivatives = true)
            : cutoff_function(std::move(cutoff_function)),
              calculate_derivatives(calculate_derivatives) {
        }
        ThreeBodyDescriptorFinder(double cutoff, double cutoff_transition_width = 0.5, bool calculate_derivatives = true)
            : calculate_derivatives(calculate_derivatives)
        {
            cutoff_function = std::make_shared<CosCutoff>(cutoff, cutoff_transition_width);
        }
        ~ThreeBodyDescriptorFinder() = default;

        CutoffRanges getCutoff() const { return { .three_body = cutoff_function->getCutoff() }; }

        virtual ThreeBodyDescriptorsFiltered findSeparations(const AtomicStructure &atomic_structure) const;
        virtual ThreeBodyDescriptorsFiltered find(const AtomicStructure &atomic_structure) const;

        virtual ThreeBodyDescriptorsFiltered toInvariantTriplets(const ThreeBodyDescriptorsFiltered& separations) const;
        virtual NewThreeBodyDescriptor toInvariantTriplet(const NewThreeBodyDescriptor& separations) const;
        virtual ThreeBodySeparationsTransformed toInvariantTriplet(const std::array<double, 3>& separations) const;

        //void tabulate(TabulationData &table) override;
    protected:
        //std::shared_ptr<CutoffFunction> optional_angular_cutoff_ = nullptr_t;
    };
}

#endif
