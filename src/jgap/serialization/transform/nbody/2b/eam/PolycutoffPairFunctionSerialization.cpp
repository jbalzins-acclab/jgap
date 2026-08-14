#include "PolycutoffPairFunctionSerialization.hpp"
#include "jgap/core/transform/nbody/2b/eam/PolycutoffPairFunction.hpp"
#include "jgap/serialization/SerializationNode.hpp"
#include "../../../../../core/io/log/CurrentLogger.hpp"

namespace jgap {

    bool PolycutoffPairFunctionSerialization::serialize(const ValuePtr<TwoBodyTransformation<1>>& obj, SerializationNode& node) const {
        if (auto derived = obj.as<PolycutoffPairFunction>()) {
            node.writeAttribute("name", "PolycutoffPairFunction");
            node.writeAttribute("cutoff", derived->getCutoff());
            node.writeAttribute("r_min", derived->getRMin());
            node.writeAttribute("prefactor", derived->getPrefactor());
            return true;
        }
        return false;
    }

    ValuePtr<TwoBodyTransformation<1>> PolycutoffPairFunctionSerialization::deserialize(const SerializationNode& node) const {
        if (node.readOptionalStringAttribute("name") != "PolycutoffPairFunction") {
            return nullptr;
        }
        auto cutoff = node.readDoubleAttribute("cutoff");
        auto r_min = node.readDoubleAttribute("r_min");
        auto prefactor = node.readDoubleAttribute("prefactor");

        return ValuePtr<TwoBodyTransformation<1>>(PolycutoffPairFunction(cutoff, r_min, prefactor));
    }

    REGISTER_SERIALIZATION(PolycutoffPairFunctionSerialization, TwoBodyTransformation<1>);
}