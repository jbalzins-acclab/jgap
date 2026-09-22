#ifndef JGAP_COSCUTOFFSERIALIZATION_HPP
#define JGAP_COSCUTOFFSERIALIZATION_HPP
#include "jgap/core/ValuePtr.hpp"
#include "jgap/core/cutoff/CutoffFunction.hpp"
#include "jgap/serialization/Serialization.hpp"

namespace jgap {
    class CosCutoffSerialization : public Serialization<CutoffFunction> {
    public:
        bool serialize(const ValuePtr<CutoffFunction>& obj, SerializationNode& node) const override;
        ValuePtr<CutoffFunction> deserialize(const SerializationNode& node) const override;
    };
}

#endif
