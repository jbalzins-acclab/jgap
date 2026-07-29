#ifndef JGAP_WENDLANDFUNCTIONSERIALIZATION_HPP
#define JGAP_WENDLANDFUNCTIONSERIALIZATION_HPP

#include "jgap/serialization/Serialization.hpp"
#include "jgap/core/cutoff/CutoffFunction.hpp"

namespace jgap {

    class WendlandFunctionSerialization : public Serialization<CutoffFunction> {
    public:
        bool serialize(const ValuePtr<CutoffFunction>& base_obj, SerializationNode& node) const override;
        ValuePtr<CutoffFunction> deserialize(const SerializationNode& node) const override;
    };
}

#endif // JGAP_WENDLANDFUNCTIONSERIALIZATION_HPP
