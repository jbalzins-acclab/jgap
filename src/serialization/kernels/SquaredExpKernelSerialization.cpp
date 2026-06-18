#include "SquaredExpKernelSerialization.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {

    template<size_t ExpDimensions, size_t CutoffDimensions>
    void SquaredExpKernelSerialization<ExpDimensions, CutoffDimensions>::serialize(
        const SquaredExpKernel<ExpDimensions, CutoffDimensions>& kernel, SerializationNode& node) {
        node.writeAttribute("name", "SquaredExpKernel");
        node.writeAttribute<int>("exp_dimensions", ExpDimensions);
        node.writeAttribute<int>("cutoff_dimensions", CutoffDimensions);
        node.writeAttribute("energy_scale", kernel.getEnergyScale());
        node.writeDataSet("length_scales", kernel.getLengthScales());
    }

    template<size_t ExpDimensions, size_t CutoffDimensions>
    SquaredExpKernel<ExpDimensions, CutoffDimensions>
    SquaredExpKernelSerialization<ExpDimensions, CutoffDimensions>::deserialize(const SerializationNode& node) {
        if (node.readOptionalAttribute<std::string>("name") != "SquaredExpKernel") {
            JGAP_LOG_AND_THROW("Node does not contain a SquaredExpKernel");
        }

        auto exp_dimensions = node.readAttribute<int>("exp_dimensions");
        auto cutoff_dimensions = node.readAttribute<int>("cutoff_dimensions");
        if (exp_dimensions != static_cast<int>(ExpDimensions) ||
            cutoff_dimensions != static_cast<int>(CutoffDimensions)) {
            JGAP_LOG_AND_THROW(
                "SquaredExpKernel dimension mismatch: stored <{}, {}> but expected <{}, {}>",
                exp_dimensions, cutoff_dimensions, ExpDimensions, CutoffDimensions);
        }

        auto energy_scale = node.readAttribute<Real>("energy_scale");
        auto length_scales = node.readDataSet<std::array<Real, ExpDimensions>>("length_scales");

        return SquaredExpKernel<ExpDimensions, CutoffDimensions>(energy_scale, length_scales);
    }

    template class SquaredExpKernelSerialization<1, 0>;
    template class SquaredExpKernelSerialization<1, 1>;
    template class SquaredExpKernelSerialization<3, 1>;
}
