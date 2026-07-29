#ifndef JGAP_COORDINATIONSTRANSFORMATIONSERIALIZATION_HPP
#define JGAP_COORDINATIONSTRANSFORMATIONSERIALIZATION_HPP

#include "jgap/serialization/Serialization.hpp"
#include "jgap/experimental/transform/nbody/2b/CoordinationTransformation.hpp"

namespace jgap {

    template<size_t Dim>
    class CoordinationTransformationSerialization : public Serialization<TwoBodyTransformation<Dim>> {
    public:
        bool serialize(const ValuePtr<TwoBodyTransformation<Dim>>& base_obj, SerializationNode& node) const override;
        ValuePtr<TwoBodyTransformation<Dim>> deserialize(const SerializationNode& node) const override;
    };
}

#endif // JGAP_COORDINATIONSTRANSFORMATIONSERIALIZATION_HPP
