#include "CosCutoffSerialization.hpp"
#include "core/cutoff/CosCutoff.hpp"
#include "io/log/CurrentLogger.hpp"
#include "serialization/SerializationNode.hpp"

namespace jgap {

    bool CosCutoffSerialization::serialize(const ValuePtr<CutoffFunction>& obj, SerializationNode& node) const {
        auto casted_fn = dynamic_cast<const CosCutoff*>(obj.get());
        if (!casted_fn) {
            return false;
        }
        node.writeAttribute("type", "CosCutoff");
        node.writeAttribute("cutoff", casted_fn->getCutoff());
        node.writeAttribute("cutoff_transition_width", casted_fn->getCutoffTransitionWidth());
        return true;
    }

    ValuePtr<CutoffFunction> CosCutoffSerialization::deserialize(const SerializationNode& node) const {
        if (node.readOptionalAttribute<std::string>("type") != "CosCutoff") {
            return nullptr;
        }

        auto cutoff = node.readAttribute<Real>("cutoff");
        auto cutoff_transition_width = node.readAttribute<Real>("cutoff_transition_width");

        return CosCutoff(cutoff, cutoff_transition_width);
    }

}
