#ifndef JGAP_TWOBODYSUMSERIALIZATION_HPP
#define JGAP_TWOBODYSUMSERIALIZATION_HPP

#include "jgap/core/ValuePtr.hpp"
#include "jgap/core/transform/manybody/NBodyAggregator.hpp"
#include "jgap/serialization/Serialization.hpp"

namespace jgap {
    template<size_t Dim>
    class TwoBodySumSerialization : public Serialization<NBodyAggregator<Dim>> {
    public:
        bool serialize(const ValuePtr<NBodyAggregator<Dim>>& obj, SerializationNode& node) const override;
        ValuePtr<NBodyAggregator<Dim>> deserialize(const SerializationNode& node) const override;
    };
}

#endif
