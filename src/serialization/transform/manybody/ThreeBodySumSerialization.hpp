#ifndef JGAP_THREEBODYSUMSERIALIZATION_HPP
#define JGAP_THREEBODYSUMSERIALIZATION_HPP

#include "core/ValuePtr.hpp"
#include "core/transform/manybody/NBodyAggregator.hpp"
#include "serialization/Serialization.hpp"

namespace jgap {
    template<size_t Dim>
    class ThreeBodySumSerialization : public Serialization<NBodyAggregator<Dim>> {
    public:
        bool serialize(const ValuePtr<NBodyAggregator<Dim>>& obj, SerializationNode& node) const override;
        ValuePtr<NBodyAggregator<Dim>> deserialize(const SerializationNode& node) const override;
    };
}

#endif
