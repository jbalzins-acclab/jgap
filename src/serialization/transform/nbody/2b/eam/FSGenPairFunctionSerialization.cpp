#include "FSGenPairFunctionSerialization.hpp"
#include "core/transform/nbody/2b/eam/FSGenPairFunction.hpp"
#include "serialization/SerializationNode.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {

    bool FSGenPairFunctionSerialization::serialize(const ValuePtr<TwoBodyTransformation<1>>& obj, SerializationNode& node) const {
        if (auto derived = obj.as<FSGenPairFunction>()) {
            node.writeAttribute("name", "FSGenPairFunction");
            node.writeAttribute("cutoff", derived->getCutoff());
            node.writeAttribute("degree", derived->getDegree());
            node.writeAttribute("prefactor", derived->getPrefactor());
            return true;
        }
        return false;
    }

    ValuePtr<TwoBodyTransformation<1>> FSGenPairFunctionSerialization::deserialize(const SerializationNode& node) const {
        if (node.readAttribute<std::string>("name") != "FSGenPairFunction") {
            return nullptr;
        }
        auto cutoff = node.readAttribute<Real>("cutoff");
        auto degree = node.readAttribute<Real>("degree");
        auto prefactor = node.readAttribute<Real>("prefactor");

        return ValuePtr<TwoBodyTransformation<1>>(FSGenPairFunction(cutoff, degree, prefactor));
    }

    REGISTER_SERIALIZATION(FSGenPairFunctionSerialization, TwoBodyTransformation<1>);
}