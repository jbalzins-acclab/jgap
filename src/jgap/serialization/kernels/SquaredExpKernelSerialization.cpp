#include "SquaredExpKernelSerialization.hpp"
#include "jgap/io/log/CurrentLogger.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"

namespace jgap {

    template<size_t ExpDimensions, size_t CutoffDimension>
    bool SquaredExpKernelSerialization<ExpDimensions, CutoffDimension>::serialize(
        const ValuePtr<Kernel<Dim>>& obj, SerializationNode& node) const {
        const auto derived = obj.template as<KernelT>();
        if (!derived) {
            return false;
        }

        node.writeAttribute("name", "SquaredExpKernel");
        node.writeAttribute("exp_dimensions", ExpDimensions);
        node.writeAttribute("cutoff_dimensions", CutoffDimension);
        node.writeAttribute("energy_scale", derived->getEnergyScale());
        node.writeDataSet("length_scales", derived->getLengthScales());
        return true;
    }

    template<size_t ExpDimensions, size_t CutoffDimension>
    auto SquaredExpKernelSerialization<ExpDimensions, CutoffDimension>::deserialize(const SerializationNode& node) const
        -> ValuePtr<Kernel<Dim>> {
        if (node.readOptionalStringAttribute("name") != "SquaredExpKernel") {
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

    template class SquaredExpKernelSerialization<1, 0>;
    template class SquaredExpKernelSerialization<1, 1>;
    template class SquaredExpKernelSerialization<2, 0>;
    template class SquaredExpKernelSerialization<2, 1>;
    template class SquaredExpKernelSerialization<3, 0>;
    template class SquaredExpKernelSerialization<3, 1>;
    template class SquaredExpKernelSerialization<4, 0>;
    template class SquaredExpKernelSerialization<4, 1>;

    using SqExpSer_1_0 = SquaredExpKernelSerialization<1, 0>;
    using SqExpSer_1_1 = SquaredExpKernelSerialization<1, 1>;
    using SqExpSer_2_0 = SquaredExpKernelSerialization<2, 0>;
    using SqExpSer_2_1 = SquaredExpKernelSerialization<2, 1>;
    using SqExpSer_3_0 = SquaredExpKernelSerialization<3, 0>;
    using SqExpSer_3_1 = SquaredExpKernelSerialization<3, 1>;
    using SqExpSer_4_0 = SquaredExpKernelSerialization<4, 0>;
    using SqExpSer_4_1 = SquaredExpKernelSerialization<4, 1>;

    REGISTER_SERIALIZATION(SqExpSer_1_0, Kernel<1>);
    REGISTER_SERIALIZATION(SqExpSer_1_1, Kernel<2>);
    REGISTER_SERIALIZATION(SqExpSer_2_0, Kernel<2>);
    REGISTER_SERIALIZATION(SqExpSer_2_1, Kernel<3>);
    REGISTER_SERIALIZATION(SqExpSer_3_0, Kernel<3>);
    REGISTER_SERIALIZATION(SqExpSer_3_1, Kernel<4>);
    REGISTER_SERIALIZATION(SqExpSer_4_0, Kernel<4>);
    REGISTER_SERIALIZATION(SqExpSer_4_1, Kernel<5>);
}
