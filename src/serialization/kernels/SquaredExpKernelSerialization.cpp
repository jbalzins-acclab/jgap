#include "SquaredExpKernelSerialization.hpp"

namespace jgap {

    template<size_t ExpDimensions, size_t CutoffDimension>
    void SquaredExpKernelSerialization<ExpDimensions, CutoffDimension>::serialize(
        const SquaredExpKernel<ExpDimensions, CutoffDimension>& kernel, SerializationNode& node) {
        node.writeAttribute("name", "SquaredExpKernel");
        node.writeAttribute<int>("exp_dimensions", ExpDimensions);
        node.writeAttribute<int>("cutoff_dimensions", CutoffDimension);
        node.writeAttribute("energy_scale", kernel.getEnergyScale());
        node.writeDataSet("length_scales", kernel.getLengthScales());
    }

    template<size_t ExpDimensions, size_t CutoffDimension>
    SquaredExpKernel<ExpDimensions, CutoffDimension>
    SquaredExpKernelSerialization<ExpDimensions, CutoffDimension>::deserialize(const SerializationNode& node) {
        if (node.readOptionalAttribute<std::string>("name") != "SquaredExpKernel") {
            JGAP_LOG_AND_THROW("Node does not contain a SquaredExpKernel");
        }

        auto exp_dimensions = node.readAttribute<int>("exp_dimensions");
        auto cutoff_dimensions = node.readAttribute<int>("cutoff_dimensions");
        if (exp_dimensions != static_cast<int>(ExpDimensions) ||
            cutoff_dimensions != static_cast<int>(CutoffDimension)) {
            JGAP_LOG_AND_THROW("SquaredExpKernel dimension mismatch: stored <{}, {}> but expected <{}, {}>",
                               exp_dimensions, cutoff_dimensions, ExpDimensions, CutoffDimension);
        }

        auto energy_scale = node.readAttribute<Real>("energy_scale");
        auto length_scales = node.readDataSet<std::array<Real, ExpDimensions>>("length_scales");

        return SquaredExpKernel<ExpDimensions, CutoffDimension>(energy_scale, length_scales);
    }

    template class SquaredExpKernelSerialization<1, 0>;
    template class SquaredExpKernelSerialization<1, 1>;
    template class SquaredExpKernelSerialization<2, 0>;
    template class SquaredExpKernelSerialization<2, 1>;
    template class SquaredExpKernelSerialization<3, 0>;
    template class SquaredExpKernelSerialization<3, 1>;
    template class SquaredExpKernelSerialization<4, 0>;
    template class SquaredExpKernelSerialization<4, 1>;
}
