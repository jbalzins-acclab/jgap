#include "Angle3bTransformationSerialization.hpp"
#include "jgap/core/cutoff/CutoffFunction.hpp"
#include "jgap/core/transform/nbody/3b/Angle3bTransformation.hpp"
#include "jgap/io/log/CurrentLogger.hpp"
#include "jgap/serialization/SerializationNode.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"

namespace jgap {

    bool Angle3bTransformationSerialization::serialize(const ValuePtr<ThreeBodyTransformation<4>>& obj,
                                                       SerializationNode& node) const {
        if (auto derived = obj.as<Angle3bTransformation>()) {
            node.writeAttribute("name", "Angle3bTransformation");
            auto cutoff_node = node.createGroup("cutoff");
            SerializationRegistry<CutoffFunction>::serialize(derived->getCutoffFunction(), cutoff_node);
            return true;
        }
        return false;
    }

    ValuePtr<ThreeBodyTransformation<4>> Angle3bTransformationSerialization::deserialize(
        const SerializationNode& node) const {
        if (node.readOptionalStringAttribute("name") != "Angle3bTransformation") {
            return nullptr;
        }
        auto cutoff_node_opt = node.getGroup("cutoff");
        if (!cutoff_node_opt) {
            JGAP_LOG_AND_THROW("Missing 'cutoff' group in Angle3bTransformation serialization");
        }
        auto cutoff = SerializationRegistry<CutoffFunction>::deserialize(cutoff_node_opt.value());
        return ValuePtr<ThreeBodyTransformation<4>>(Angle3bTransformation(cutoff));
    }

    REGISTER_SERIALIZATION(Angle3bTransformationSerialization, ThreeBodyTransformation<4>);
}
