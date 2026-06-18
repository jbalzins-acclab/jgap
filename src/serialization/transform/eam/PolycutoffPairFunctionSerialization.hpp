#ifndef JGAP_POLYCUTOFFPAIRFUNCTIONSERIALIZATION_HPP
#define JGAP_POLYCUTOFFPAIRFUNCTIONSERIALIZATION_HPP

#include "core/transform/ClusterTransformation.hpp"
#include "serialization/Serialization.hpp"
#include "serialization/SerializationRegistry.hpp"
#include "core/ValuePtr.hpp"

namespace jgap {
    class PolycutoffPairFunctionSerialization : public Serialization<ClusterTransformation<1, 2>> {
    public:
        bool serialize(const ValuePtr<ClusterTransformation<1, 2>>& obj, SerializationNode& node) const override;
        ValuePtr<ClusterTransformation<1, 2>> deserialize(const SerializationNode& node) const override;
    };
}

#endif