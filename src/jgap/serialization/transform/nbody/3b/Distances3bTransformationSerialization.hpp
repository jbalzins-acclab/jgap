#ifndef JGAP_DISTANCES3BTRANSFORMATIONSERIALIZATION_HPP
#define JGAP_DISTANCES3BTRANSFORMATIONSERIALIZATION_HPP

#include "jgap/core/ValuePtr.hpp"
#include "jgap/core/transform/nbody/3b/ThreeBodyTransformation.hpp"
#include "jgap/serialization/Serialization.hpp"

namespace jgap {
    class Distances3bTransformationSerialization : public Serialization<ThreeBodyTransformation<4>> {
    public:
        bool serialize(const ValuePtr<ThreeBodyTransformation<4>>& obj, SerializationNode& node) const override;
        ValuePtr<ThreeBodyTransformation<4>> deserialize(const SerializationNode& node) const override;
    };
}

#endif
