#include "MeamTransformationSerialization.hpp"
#include "jgap/core/cutoff/CutoffFunction.hpp"
#include "jgap/io/log/CurrentLogger.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"

namespace jgap {

    bool MeamTransformationSerialization::serialize(const ValuePtr<ThreeBodyTransformation<3>>& obj,
                                                    SerializationNode& node) const {
        if (auto derived = obj.as<MeamTransformation>()) {
            node.writeAttribute("name", "MeamTransformation");
            auto cutoff_node = node.createGroup("cutoff");
            SerializationRegistry<CutoffFunction>::serialize(derived->getCutoffFunction(), cutoff_node);
            return true;
        }
        return false;
    }

    ValuePtr<ThreeBodyTransformation<3>> MeamTransformationSerialization::deserialize(
        const SerializationNode& node) const {
        if (node.readOptionalStringAttribute("name") != "MeamTransformation") {
            return nullptr;
        }
        auto cutoff_node_opt = node.getGroup("cutoff");
        if (!cutoff_node_opt) {
            JGAP_LOG_AND_THROW("Missing 'cutoff' group in MeamTransformation serialization");
        }
        auto cutoff = SerializationRegistry<CutoffFunction>::deserialize(cutoff_node_opt.value());
        return ValuePtr<ThreeBodyTransformation<3>>(MeamTransformation(cutoff));
    }

    REGISTER_SERIALIZATION(MeamTransformationSerialization, ThreeBodyTransformation<3>);
}
