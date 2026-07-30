#ifndef JGAP_MEAMTRANSFORMATIONSERIALIZATION_HPP
#define JGAP_MEAMTRANSFORMATIONSERIALIZATION_HPP

#include "jgap/experimental/transform/nbody/3b/MeamTransformation.hpp"
#include "jgap/serialization/Serialization.hpp"
#include "jgap/serialization/SerializationNode.hpp"

namespace jgap {

    class MeamTransformationSerialization final : public Serialization<ThreeBodyTransformation<3>> {
    public:
        bool serialize(const ValuePtr<ThreeBodyTransformation<3>>& obj, SerializationNode& node) const override;
        ValuePtr<ThreeBodyTransformation<3>> deserialize(const SerializationNode& node) const override;
    };

}

#endif // JGAP_MEAMTRANSFORMATIONSERIALIZATION_HPP
