#ifndef JGAP_ANGLE3BTRANSFORMATIONSERIALIZATION_HPP
#define JGAP_ANGLE3BTRANSFORMATIONSERIALIZATION_HPP
#include "core/transform/NBodyTransformation.hpp"
#include "serialization/Serialization.hpp"
#include "serialization/SerializationRegistry.hpp"
#include "core/ValuePtr.hpp"

namespace jgap {
    class Angle3bTransformationSerialization : public Serialization<NBodyTransformation<4, 3>> {
    public:
        bool serialize(const ValuePtr<NBodyTransformation<4, 3>>& obj, SerializationNode& node) const override;
        ValuePtr<NBodyTransformation<4, 3>> deserialize(const SerializationNode& node) const override;
    };

}

#endif