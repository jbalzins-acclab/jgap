#include "PairDistanceTransformationSerialization.hpp"
#include "jgap/core/cutoff/CutoffFunction.hpp"
#include "jgap/core/transform/nbody/2b/PairDistanceTransformation.hpp"
#include "jgap/core/transform/nbody/2b/TwoBodyTransformation.hpp"
#include "../../../../core/io/log/CurrentLogger.hpp"
#include "jgap/serialization/SerializationNode.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"

namespace jgap {

    bool PairDistanceTransformationSerialization::serialize(const ValuePtr<TwoBodyTransformation<2>>& obj,
                                                            SerializationNode& node) const {
        if (auto derived = obj.as<PairDistanceTransformation>()) {
            node.writeAttribute("name", "TwoBodyTransformation");
            auto cutoff_node = node.createGroup("cutoff");
            SerializationRegistry<CutoffFunction>::serialize(derived->getCutoffFunction(), cutoff_node);
            return true;
        }
        return false;
    }

    ValuePtr<TwoBodyTransformation<2>> PairDistanceTransformationSerialization::deserialize(
        const SerializationNode& node) const {
        if (node.readOptionalStringAttribute("name") != "TwoBodyTransformation") {
            return nullptr;
        }
        auto cutoff_node_opt = node.getGroup("cutoff");
        if (!cutoff_node_opt) {
            JGAP_LOG_AND_THROW("Missing 'cutoff' group in TwoBodyTransformation serialization");
        }
        auto cutoff = SerializationRegistry<CutoffFunction>::deserialize(cutoff_node_opt.value());

        return ValuePtr<TwoBodyTransformation<2>>(PairDistanceTransformation(cutoff));
    }

    REGISTER_SERIALIZATION(PairDistanceTransformationSerialization, TwoBodyTransformation<2>);
}
