#ifndef JGAP_TWOBODYTRANSFORMATIONSERIALIZATION_HPP
#define JGAP_TWOBODYTRANSFORMATIONSERIALIZATION_HPP
#include "core/transform/NBodyTransformation.hpp"
#include "serialization/Serialization.hpp"
#include "serialization/SerializationRegistry.hpp"
#include "core/ValuePtr.hpp"

namespace jgap {
    class TwoBodyTransformationSerialization : public Serialization<NBodyTransformation<2, 2>> {
    public:
        bool serialize(const ValuePtr<NBodyTransformation<2, 2>>& obj, SerializationNode& node) const override;
        ValuePtr<NBodyTransformation<2, 2>> deserialize(const SerializationNode& node) const override;
    };

}

#endif