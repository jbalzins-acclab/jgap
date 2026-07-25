#ifndef JGAP_WENDLANDKERNELSERIALIZATION_HPP
#define JGAP_WENDLANDKERNELSERIALIZATION_HPP

#include "jgap/experimental/kernels/WendlandKernel.hpp"
#include "jgap/serialization/Serialization.hpp"
#include "jgap/serialization/SerializationNode.hpp"

namespace jgap {

    template<size_t ExpDimensions, size_t CutoffDimension>
    class WendlandKernelSerialization : public Serialization<Kernel<ExpDimensions + CutoffDimension>> {
    public:
        static constexpr size_t Dim = ExpDimensions + CutoffDimension;
        using KernelT = WendlandKernel<ExpDimensions, CutoffDimension>;

        bool serialize(const ValuePtr<Kernel<Dim>>& obj, SerializationNode& node) const override;
        ValuePtr<Kernel<Dim>> deserialize(const SerializationNode& node) const override;
    };
}

#endif
