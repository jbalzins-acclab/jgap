#ifndef JGAP_SPLINEPAIRTRANSFORMATIONSERIALIZATION_HPP
#define JGAP_SPLINEPAIRTRANSFORMATIONSERIALIZATION_HPP

#include "jgap/core/transform/nbody/2b/TwoBodyTransformation.hpp"
#include "jgap/serialization/Serialization.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"
#include "jgap/core/ValuePtr.hpp"

namespace jgap {
    class SplinePairTransformationSerialization : public Serialization<TwoBodyTransformation<1>> {
    public:
        bool serialize(const ValuePtr<TwoBodyTransformation<1>>& obj, SerializationNode& node) const override;
        ValuePtr<TwoBodyTransformation<1>> deserialize(const SerializationNode& node) const override;
    };
}

#endif