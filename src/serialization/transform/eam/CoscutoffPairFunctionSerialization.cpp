#include "CoscutoffPairFunctionSerialization.hpp"
#include "core/transform/eam/CoscutoffPairFunction.hpp"
#include "serialization/SerializationNode.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {

    bool CoscutoffPairFunctionSerialization::serialize(const ValuePtr<ClusterTransformation<1, 2>>& obj, SerializationNode& node) const {
        if (auto derived = obj.as<CoscutoffPairFunction>()) {
            node.writeAttribute("name", "CoscutoffPairFunction");
            node.writeAttribute("cutoff", derived->getCutoff());
            node.writeAttribute("r_min", derived->getRMin());
            node.writeAttribute("prefactor", derived->getPrefactor());
            return true;
        }
        return false;
    }

    ValuePtr<ClusterTransformation<1, 2>> CoscutoffPairFunctionSerialization::deserialize(const SerializationNode& node) const {
        if (node.readAttribute<std::string>("name") != "CoscutoffPairFunction") {
            return nullptr;
        }
        auto cutoff = node.readAttribute<Real>("cutoff");
        auto r_min = node.readAttribute<Real>("r_min");
        auto prefactor = node.readAttribute<Real>("prefactor");

        return ValuePtr<ClusterTransformation<1, 2>>(CoscutoffPairFunction(cutoff, r_min, prefactor));
    }

    REGISTER_SERIALIZATION(CoscutoffPairFunctionSerialization, ClusterTransformation<1, 2>);
}