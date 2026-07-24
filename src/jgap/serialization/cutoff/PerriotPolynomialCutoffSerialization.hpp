#ifndef JGAP_PERRIOTPOLYNOMIALCUTOFFSERIALIZATION_HPP
#define JGAP_PERRIOTPOLYNOMIALCUTOFFSERIALIZATION_HPP
#include "jgap/core/ValuePtr.hpp"
#include "jgap/core/cutoff/CutoffFunction.hpp"
#include "jgap/serialization/Serialization.hpp"

namespace jgap {
    class PerriotPolynomialCutoffSerialization : public Serialization<CutoffFunction> {
    public:
        bool serialize(const ValuePtr<CutoffFunction>& obj, SerializationNode& node) const override;
        ValuePtr<CutoffFunction> deserialize(const SerializationNode& node) const override;
    };
}

#endif
