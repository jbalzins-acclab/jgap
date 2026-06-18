#include "TwoBodyTransformationSerialization.hpp"
#include "core/transform/2b/TwoBodyTransformation.hpp"
#include "serialization/SerializationNode.hpp"
#include "core/cutoff/CutoffFunction.hpp"
#include "serialization/SerializationRegistry.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {

    bool TwoBodyTransformationSerialization::serialize(const ValuePtr<ClusterTransformation<2, 2>>& obj,
        SerializationNode& node) const
    {
        if (auto derived = obj.as<TwoBodyTransformation>()) {
            node.writeAttribute("name", "TwoBodyTransformation");
            auto cutoff_node = node.createGroup("cutoff");
            SerializationRegistry<CutoffFunction>::serialize(derived->getCutoff(), cutoff_node);
            return true;
        }
        return false;
    }

    ValuePtr<ClusterTransformation<2, 2>> TwoBodyTransformationSerialization::deserialize(
        const SerializationNode& node) const
    {
        if (node.readOptionalAttribute<std::string>("name") != "TwoBodyTransformation") {
            return nullptr;
        }
        auto cutoff_node_opt = node.getGroup("cutoff");
        if (!cutoff_node_opt) {
            JGAP_LOG_AND_THROW("Missing 'cutoff' group in TwoBodyTransformation serialization");
        }
        auto cutoff = SerializationRegistry<CutoffFunction>::deserialize(cutoff_node_opt.value());
        return ValuePtr<ClusterTransformation<2, 2>>(TwoBodyTransformation(cutoff));
    }

    REGISTER_SERIALIZATION(TwoBodyTransformationSerialization, ClusterTransformation<2, 2>);
}