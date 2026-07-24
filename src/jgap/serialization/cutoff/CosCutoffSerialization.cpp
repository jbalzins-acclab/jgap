#include "CosCutoffSerialization.hpp"
#include "jgap/core/cutoff/CosCutoff.hpp"
#include "jgap/serialization/SerializationNode.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"

namespace jgap {

    bool CosCutoffSerialization::serialize(const ValuePtr<CutoffFunction>& obj, SerializationNode& node) const {
        auto casted_fn = obj.as<CosCutoff>();
        if (!casted_fn) {
            return false;
        }
        node.writeAttribute("type", "CosCutoff");
        node.writeAttribute("cutoff", casted_fn->getCutoff());
        node.writeAttribute("cutoff_transition_width", casted_fn->getCutoffTransitionWidth());
        return true;
    }

    ValuePtr<CutoffFunction> CosCutoffSerialization::deserialize(const SerializationNode& node) const {
        if (node.readOptionalStringAttribute("type") != "CosCutoff") {
            return nullptr;
        }

        auto cutoff = node.readDoubleAttribute("cutoff");
        auto cutoff_transition_width = node.readDoubleAttribute("cutoff_transition_width");

        return CosCutoff(cutoff, cutoff_transition_width);
    }

    REGISTER_SERIALIZATION(CosCutoffSerialization, CutoffFunction)
}
