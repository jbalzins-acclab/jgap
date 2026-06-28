#ifndef JGAP_SPLINEPAIRTRANSFORMATIONSERIALIZATION_HPP
#define JGAP_SPLINEPAIRTRANSFORMATIONSERIALIZATION_HPP

#include "core/transform/NBodyTransformation.hpp"
#include "serialization/Serialization.hpp"
#include "serialization/SerializationRegistry.hpp"
#include "core/ValuePtr.hpp"

namespace jgap {
    class SplinePairTransformationSerialization : public Serialization<NBodyTransformation<1, 2>> {
    public:
        bool serialize(const ValuePtr<NBodyTransformation<1, 2>>& obj, SerializationNode& node) const override;
        ValuePtr<NBodyTransformation<1, 2>> deserialize(const SerializationNode& node) const override;
    };
}

#endif