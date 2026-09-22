#ifndef JGAP_CAUCHYKERNELSERIALIZATION_HPP
#define JGAP_CAUCHYKERNELSERIALIZATION_HPP

#include "jgap/experimental/kernels/CauchyKernel.hpp"
#include "jgap/serialization/Serialization.hpp"
#include "jgap/serialization/SerializationNode.hpp"

namespace jgap {

    template<size_t ExpDimensions, size_t CutoffDimension>
    class CauchyKernelSerialization : public Serialization<Kernel<ExpDimensions + CutoffDimension>> {
    public:
        static constexpr size_t Dim = ExpDimensions + CutoffDimension;
        using KernelT = CauchyKernel<ExpDimensions, CutoffDimension>;

        bool serialize(const ValuePtr<Kernel<Dim>>& obj, SerializationNode& node) const override;
        ValuePtr<Kernel<Dim>> deserialize(const SerializationNode& node) const override;
    };
}

#endif
