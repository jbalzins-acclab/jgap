#include "PerriotPolynomialCutoffSerialization.hpp"
#include "core/cutoff/PerriotPolynomialCutoff.hpp"
#include "io/log/CurrentLogger.hpp"
#include "serialization/SerializationNode.hpp"

namespace jgap {

    bool PerriotPolynomialCutoffSerialization::serialize(const ValuePtr<CutoffFunction>& obj,
                                                         SerializationNode& node) const {
        auto casted_fn = dynamic_cast<const PerriotPolynomialCutoff*>(obj.get());
        if (!casted_fn) {
            return false;
        }
        node.writeAttribute("type", "PerriotPolynomialCutoff");
        node.writeAttribute("cutoff", casted_fn->getCutoff());
        node.writeAttribute("cutoff_transition_width", casted_fn->getCutoffTransitionWidth());
        return true;
    }

    ValuePtr<CutoffFunction> PerriotPolynomialCutoffSerialization::deserialize(const SerializationNode& node) const {
        if (node.readOptionalAttribute<std::string>("type") != "PerriotPolynomialCutoff") {
            return nullptr;
        }

        auto cutoff = node.readAttribute<Real>("cutoff");
        auto cutoff_transition_width = node.readAttribute<Real>("cutoff_transition_width");

        return PerriotPolynomialCutoff(cutoff, cutoff_transition_width);
    }

}
