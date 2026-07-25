#include "CutoffJK3bTransformationSerialization.hpp"

#include "jgap/core/cutoff/CutoffFunction.hpp"
#include "jgap/experimental/transform/nbody/3b/CutoffJK3bTransformation.hpp"
#include "jgap/io/log/CurrentLogger.hpp"
#include "jgap/serialization/SerializationNode.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"

namespace jgap {

    bool CutoffJK3bTransformationSerialization::serialize(const ValuePtr<ThreeBodyTransformation<4>>& obj,
                                                           SerializationNode& node) const {
        if (auto derived = obj.as<CutoffJK3bTransformation>()) {
            node.writeAttribute("name", "CutoffJK3bTransformation");
            auto main_cutoff_node = node.createGroup("main_cutoff");
            SerializationRegistry<CutoffFunction>::serialize(derived->getMainCutoffFunction(), main_cutoff_node);
            auto cutoff_12_node = node.createGroup("cutoff_12");
            SerializationRegistry<CutoffFunction>::serialize(derived->getCutoff12Function(), cutoff_12_node);
            return true;
        }
        return false;
    }

    ValuePtr<ThreeBodyTransformation<4>> CutoffJK3bTransformationSerialization::deserialize(
        const SerializationNode& node) const {
        if (node.readOptionalStringAttribute("name") != "CutoffJK3bTransformation") {
            return nullptr;
        }
        auto main_cutoff_node_opt = node.getGroup("main_cutoff");
        if (!main_cutoff_node_opt) {
            JGAP_LOG_AND_THROW("Missing 'main_cutoff' group in CutoffJK3bTransformation serialization");
        }
        auto main_cutoff = SerializationRegistry<CutoffFunction>::deserialize(main_cutoff_node_opt.value());

        auto cutoff_12_node_opt = node.getGroup("cutoff_12");
        if (!cutoff_12_node_opt) {
            JGAP_LOG_AND_THROW("Missing 'cutoff_12' group in CutoffJK3bTransformation serialization");
        }
        auto cutoff_12 = SerializationRegistry<CutoffFunction>::deserialize(cutoff_12_node_opt.value());

        return ValuePtr<ThreeBodyTransformation<4>>(CutoffJK3bTransformation(main_cutoff, cutoff_12));
    }

    REGISTER_SERIALIZATION(CutoffJK3bTransformationSerialization, ThreeBodyTransformation<4>);
}
