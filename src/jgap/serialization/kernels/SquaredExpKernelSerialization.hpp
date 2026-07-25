#ifndef JGAP_SQUAREDEXPKERNELSERIALIZATION_HPP
#define JGAP_SQUAREDEXPKERNELSERIALIZATION_HPP

#include "jgap/core/kernels/SquaredExpKernel.hpp"
#include "jgap/serialization/Serialization.hpp"
#include "jgap/serialization/SerializationNode.hpp"

namespace jgap {

    template<size_t ExpDimensions, size_t CutoffDimension>
    class SquaredExpKernelSerialization : public Serialization<Kernel<ExpDimensions + CutoffDimension>> {
    public:
        static constexpr size_t Dim = ExpDimensions + CutoffDimension;
        using KernelT = SquaredExpKernel<ExpDimensions, CutoffDimension>;

        bool serialize(const ValuePtr<Kernel<Dim>>& obj, SerializationNode& node) const override;
        ValuePtr<Kernel<Dim>> deserialize(const SerializationNode& node) const override;
    };
}

#endif
