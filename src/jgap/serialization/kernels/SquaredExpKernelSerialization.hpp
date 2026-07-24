#ifndef JGAP_SQUAREDEXPKERNELSERIALIZATION_HPP
#define JGAP_SQUAREDEXPKERNELSERIALIZATION_HPP

#include "jgap/core/kernels/SquaredExpKernel.hpp"
#include "jgap/serialization/SerializationNode.hpp"

namespace jgap {

    // SquaredExpKernel is used as a value type parameterised by its dimensions rather than through a
    // Cloneable polymorphic base, so it is serialized directly instead of via SerializationRegistry.
    // Explicit instantiations live in the .cpp; add one there for every dimension combination used.
    template<size_t ExpDimensions, size_t CutoffDimension>
    class SquaredExpKernelSerialization {
    public:
        static void serialize(const SquaredExpKernel<ExpDimensions, CutoffDimension>& kernel, SerializationNode& node);
        static SquaredExpKernel<ExpDimensions, CutoffDimension> deserialize(const SerializationNode& node);
    };
}

#endif
