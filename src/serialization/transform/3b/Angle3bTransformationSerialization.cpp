#include "Angle3bTransformationSerialization.hpp"
#include "core/transform/3b/Angle3bTransformation.hpp"
#include "serialization/SerializationNode.hpp"
#include "core/cutoff/CutoffFunction.hpp"
#include "serialization/SerializationRegistry.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {

    bool Angle3bTransformationSerialization::serialize(const ValuePtr<ClusterTransformation<4, 3>>& obj,
        SerializationNode& node) const
    {
        if (auto derived = obj.as<Angle3bTransformation>()) {
            node.writeAttribute("name", "Angle3bTransformation");
            auto cutoff_node = node.createGroup("cutoff");
            SerializationRegistry<CutoffFunction>::serialize(derived->getCutoff(), cutoff_node);
            return true;
        }
        return false;
    }

    ValuePtr<ClusterTransformation<4, 3>> Angle3bTransformationSerialization::deserialize(
        const SerializationNode& node) const
    {
        if (node.readOptionalAttribute<std::string>("name") != "Angle3bTransformation") {
            return nullptr;
        }
        auto cutoff_node_opt = node.getGroup("cutoff");
        if (!cutoff_node_opt) {
            JGAP_LOG_AND_THROW("Missing 'cutoff' group in Angle3bTransformation serialization");
        }
        auto cutoff = SerializationRegistry<CutoffFunction>::deserialize(cutoff_node_opt.value());
        return ValuePtr<ClusterTransformation<4, 3>>(Angle3bTransformation(cutoff));
    }

    REGISTER_SERIALIZATION(Angle3bTransformationSerialization, ClusterTransformation<4, 3>);
}