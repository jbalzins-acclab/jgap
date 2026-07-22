#ifndef JGAP_PAIRDISTANCETRANSFORMATIONSERIALIZATION_HPP
#define JGAP_PAIRDISTANCETRANSFORMATIONSERIALIZATION_HPP

#include "core/ValuePtr.hpp"
#include "core/transform/nbody/2b/TwoBodyTransformation.hpp"
#include "serialization/Serialization.hpp"
#include "serialization/SerializationRegistry.hpp"

namespace jgap {
    class PairDistanceTransformationSerialization : public Serialization<TwoBodyTransformation<2>> {
    public:
        bool serialize(const ValuePtr<TwoBodyTransformation<2>>& obj, SerializationNode& node) const override;
        ValuePtr<TwoBodyTransformation<2>> deserialize(const SerializationNode& node) const override;
    };
}

#endif
