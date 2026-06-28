#include "PolycutoffPairFunctionSerialization.hpp"
#include "core/transform/eam/PolycutoffPairFunction.hpp"
#include "serialization/SerializationNode.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {

    bool PolycutoffPairFunctionSerialization::serialize(const ValuePtr<NBodyTransformation<1, 2>>& obj, SerializationNode& node) const {
        if (auto derived = obj.as<PolycutoffPairFunction>()) {
            node.writeAttribute("name", "PolycutoffPairFunction");
            node.writeAttribute("cutoff", derived->getCutoff());
            node.writeAttribute("r_min", derived->getRMin());
            node.writeAttribute("prefactor", derived->getPrefactor());
            return true;
        }
        return false;
    }

    ValuePtr<NBodyTransformation<1, 2>> PolycutoffPairFunctionSerialization::deserialize(const SerializationNode& node) const {
        if (node.readOptionalAttribute<std::string>("name") != "PolycutoffPairFunction") {
            return nullptr;
        }
        auto cutoff = node.readAttribute<Real>("cutoff");
        auto r_min = node.readAttribute<Real>("r_min");
        auto prefactor = node.readAttribute<Real>("prefactor");

        return ValuePtr<NBodyTransformation<1, 2>>(PolycutoffPairFunction(cutoff, r_min, prefactor));
    }

    REGISTER_SERIALIZATION(PolycutoffPairFunctionSerialization, NBodyTransformation<1, 2>);
}