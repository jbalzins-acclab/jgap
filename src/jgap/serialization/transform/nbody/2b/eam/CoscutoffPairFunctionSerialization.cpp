#include "CoscutoffPairFunctionSerialization.hpp"
#include "jgap/core/transform/nbody/2b/eam/CoscutoffPairFunction.hpp"
#include "jgap/serialization/SerializationNode.hpp"
#include "../../../../../core/io/log/CurrentLogger.hpp"

namespace jgap {

    bool CoscutoffPairFunctionSerialization::serialize(const ValuePtr<TwoBodyTransformation<1>>& obj, SerializationNode& node) const {
        if (auto derived = obj.as<CoscutoffPairFunction>()) {
            node.writeAttribute("name", "CoscutoffPairFunction");
            node.writeAttribute("cutoff", derived->getCutoff());
            node.writeAttribute("r_min", derived->getRMin());
            node.writeAttribute("prefactor", derived->getPrefactor());
            return true;
        }
        return false;
    }

    ValuePtr<TwoBodyTransformation<1>> CoscutoffPairFunctionSerialization::deserialize(const SerializationNode& node) const {
        if (node.readStringAttribute("name") != "CoscutoffPairFunction") {
            return nullptr;
        }
        auto cutoff = node.readDoubleAttribute("cutoff");
        auto r_min = node.readDoubleAttribute("r_min");
        auto prefactor = node.readDoubleAttribute("prefactor");

        return ValuePtr<TwoBodyTransformation<1>>(CoscutoffPairFunction(cutoff, r_min, prefactor));
    }

    REGISTER_SERIALIZATION(CoscutoffPairFunctionSerialization, TwoBodyTransformation<1>);
}