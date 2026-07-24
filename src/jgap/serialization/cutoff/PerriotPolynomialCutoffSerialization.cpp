#include "PerriotPolynomialCutoffSerialization.hpp"
#include "jgap/core/cutoff/PerriotPolynomialCutoff.hpp"
#include "jgap/serialization/SerializationNode.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"

namespace jgap {

    bool PerriotPolynomialCutoffSerialization::serialize(const ValuePtr<CutoffFunction>& obj,
                                                         SerializationNode& node) const {
        auto casted_fn = obj.as<PerriotPolynomialCutoff>();
        if (!casted_fn) {
            return false;
        }
        node.writeAttribute("type", "PerriotPolynomialCutoff");
        node.writeAttribute("cutoff", casted_fn->getCutoff());
        node.writeAttribute("cutoff_transition_width", casted_fn->getCutoffTransitionWidth());
        return true;
    }

    ValuePtr<CutoffFunction> PerriotPolynomialCutoffSerialization::deserialize(const SerializationNode& node) const {
        if (node.readOptionalStringAttribute("type") != "PerriotPolynomialCutoff") {
            return nullptr;
        }

        auto cutoff = node.readDoubleAttribute("cutoff");
        auto cutoff_transition_width = node.readDoubleAttribute("cutoff_transition_width");

        return PerriotPolynomialCutoff(cutoff, cutoff_transition_width);
    }

    REGISTER_SERIALIZATION(PerriotPolynomialCutoffSerialization, CutoffFunction)
}
