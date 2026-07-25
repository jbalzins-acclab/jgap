#include "WendlandKernelSerialization.hpp"
#include "jgap/io/log/CurrentLogger.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"

namespace jgap {

    template<size_t ExpDimensions, size_t CutoffDimension>
    bool WendlandKernelSerialization<ExpDimensions, CutoffDimension>::serialize(
        const ValuePtr<Kernel<Dim>>& obj, SerializationNode& node
    ) const {
        const auto derived = obj.template as<KernelT>();
        if (!derived) {
            return false;
        }

        node.writeAttribute("name", "WendlandKernel");
        node.writeAttribute("exp_dimensions", ExpDimensions);
        node.writeAttribute("cutoff_dimensions", CutoffDimension);
        node.writeAttribute("energy_scale", derived->getEnergyScale());
        node.writeDataSet("length_scales", derived->getLengthScales());
        return true;
    }

    template<size_t ExpDimensions, size_t CutoffDimension>
    auto WendlandKernelSerialization<ExpDimensions, CutoffDimension>::deserialize(const SerializationNode& node) const
        -> ValuePtr<Kernel<Dim>> {
        if (node.readOptionalStringAttribute("name") != "WendlandKernel") {
            return nullptr;
        }

        auto exp_dimensions = node.readIntAttribute("exp_dimensions");
        auto cutoff_dimensions = node.readIntAttribute("cutoff_dimensions");
        if (exp_dimensions != static_cast<int>(ExpDimensions) ||
            cutoff_dimensions != static_cast<int>(CutoffDimension)) {
            return nullptr;
        }

        auto energy_scale = node.readDoubleAttribute("energy_scale");
        auto length_scales = node.readDoubleArrayDataSet<ExpDimensions>("length_scales");

        return ValuePtr<Kernel<Dim>>(KernelT(energy_scale, length_scales));
    }

    template class WendlandKernelSerialization<1, 0>;
    template class WendlandKernelSerialization<1, 1>;
    template class WendlandKernelSerialization<2, 0>;
    template class WendlandKernelSerialization<2, 1>;
    template class WendlandKernelSerialization<3, 0>;
    template class WendlandKernelSerialization<3, 1>;
    template class WendlandKernelSerialization<4, 0>;
    template class WendlandKernelSerialization<4, 1>;

    using WendlandSer_1_0 = WendlandKernelSerialization<1, 0>;
    using WendlandSer_1_1 = WendlandKernelSerialization<1, 1>;
    using WendlandSer_2_0 = WendlandKernelSerialization<2, 0>;
    using WendlandSer_2_1 = WendlandKernelSerialization<2, 1>;
    using WendlandSer_3_0 = WendlandKernelSerialization<3, 0>;
    using WendlandSer_3_1 = WendlandKernelSerialization<3, 1>;
    using WendlandSer_4_0 = WendlandKernelSerialization<4, 0>;
    using WendlandSer_4_1 = WendlandKernelSerialization<4, 1>;

    REGISTER_SERIALIZATION(WendlandSer_1_0, Kernel<1>);
    REGISTER_SERIALIZATION(WendlandSer_1_1, Kernel<2>);
    REGISTER_SERIALIZATION(WendlandSer_2_0, Kernel<2>);
    REGISTER_SERIALIZATION(WendlandSer_2_1, Kernel<3>);
    REGISTER_SERIALIZATION(WendlandSer_3_0, Kernel<3>);
    REGISTER_SERIALIZATION(WendlandSer_3_1, Kernel<4>);
    REGISTER_SERIALIZATION(WendlandSer_4_0, Kernel<4>);
    REGISTER_SERIALIZATION(WendlandSer_4_1, Kernel<5>);
}
