#include "WendlandFunctionSerialization.hpp"
#include "jgap/experimental/cutoff/WendlandFunction.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"

namespace jgap {

    bool WendlandFunctionSerialization::serialize(
        const ValuePtr<CutoffFunction>& base_obj, SerializationNode& node
    ) const {
        if (auto obj = dynamic_cast<const WendlandFunction*>(base_obj.get())) {
            node.writeAttribute("name", "WendlandFunction");
            node.writeAttribute("r_min", obj->getRMin());
            node.writeAttribute("r_max", obj->getRMax());
            return true;
        }
        return false;
    }

    ValuePtr<CutoffFunction> WendlandFunctionSerialization::deserialize(const SerializationNode& node) const {
        if (node.readOptionalStringAttribute("name") != "WendlandFunction") {
            return nullptr;
        }

        Real r_min = node.readDoubleAttribute("r_min");
        Real r_max = node.readDoubleAttribute("r_max");

        return WendlandFunction(r_min, r_max);
    }

    REGISTER_SERIALIZATION(WendlandFunctionSerialization, CutoffFunction);
}
