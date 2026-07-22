#ifndef JGAP_ANGLE3BTRANSFORMATIONSERIALIZATION_HPP
#define JGAP_ANGLE3BTRANSFORMATIONSERIALIZATION_HPP
#include "core/transform/nbody/3b/ThreeBodyTransformation.hpp"
#include "serialization/Serialization.hpp"
#include "serialization/SerializationRegistry.hpp"
#include "core/ValuePtr.hpp"

namespace jgap {
    class Angle3bTransformationSerialization : public Serialization<ThreeBodyTransformation<4>> {
    public:
        bool serialize(const ValuePtr<ThreeBodyTransformation<4>>& obj, SerializationNode& node) const override;
        ValuePtr<ThreeBodyTransformation<4>> deserialize(const SerializationNode& node) const override;
    };

}

#endif