#ifndef JGAP_PERRIOTPOLYNOMIALCUTOFFSERIALIZATION_HPP
#define JGAP_PERRIOTPOLYNOMIALCUTOFFSERIALIZATION_HPP
#include "core/cutoff/CutoffFunction.hpp"
#include "serialization/Serialization.hpp"
#include "serialization/SerializationRegistry.hpp"
#include "core/ValuePtr.hpp"

namespace jgap {
    class PerriotPolynomialCutoffSerialization : public Serialization<CutoffFunction> {
    public:
        bool serialize(const ValuePtr<CutoffFunction>& obj, SerializationNode& node) const override;
        ValuePtr<CutoffFunction> deserialize(const SerializationNode& node) const override;
    };
}

#endif