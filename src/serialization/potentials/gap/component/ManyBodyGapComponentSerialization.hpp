#ifndef JGAP_MANYBODYGAPCOMPONENTSERIALIZATION_HPP
#define JGAP_MANYBODYGAPCOMPONENTSERIALIZATION_HPP

#include "core/potentials/gap/component/GapComponent.hpp"
#include "core/potentials/gap/component/ManyBodyGapComponent.hpp"
#include "core/transform/aggregated/TransformationAggregator.hpp"
#include "core/kernels/SquaredExpKernel.hpp"
#include "serialization/Serialization.hpp"
#include "serialization/SerializationRegistry.hpp"
#include "serialization/SerializationNode.hpp"
#include "serialization/kernels/SquaredExpKernelSerialization.hpp"
#include "core/ValuePtr.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {

    // Serializer for the ManyBodyGapComponent<Dim, TKernel> family. As with NBodyGapComponent the kernel
    // is assumed to be a SquaredExpKernel. Explicit instantiations and registrations live in the .cpp.
    template<size_t Dim, typename TKernel>
    class ManyBodyGapComponentSerialization : public Serialization<GapComponent> {
        using ComponentT = ManyBodyGapComponent<Dim, TKernel>;
        using KernelSerialization = SquaredExpKernelSerialization<TKernel::ExpDim, TKernel::CutoffDim>;

    public:
        bool serialize(const ValuePtr<GapComponent>& obj, SerializationNode& node) const override {
            auto derived = obj.template as<ComponentT>();
            if (!derived) {
                return false;
            }

            node.writeAttribute("name", "ManyBodyGapComponent");
            node.writeAttribute("dim", Dim);

            auto kernel_group = node.createGroup("kernel");
            KernelSerialization::serialize(derived->getKernel(), kernel_group);

            auto aggregator_group = node.createGroup("aggregator");
            SerializationRegistry<TransformationAggregator<Dim>>::serialize(
                derived->getAggregator(), aggregator_group);

            node.saveSparsePoints<Dim>(derived->getSparsePoints());

            const auto& coefficients = derived->getCoefficients();
            if (!coefficients.empty()) {
                node.writeDataSet("coefficients", coefficients);
            }

            return true;
        }

        ValuePtr<GapComponent> deserialize(const SerializationNode& node) const override {
            if (node.readOptionalAttribute<std::string>("name") != "ManyBodyGapComponent" ||
                node.readOptionalAttribute<size_t>("dim") != Dim) {
                return nullptr;
            }

            auto kernel_group_opt = node.getGroup("kernel");
            if (!kernel_group_opt) JGAP_LOG_AND_THROW("Missing 'kernel' group in ManyBodyGapComponent serialization");
            TKernel kernel = KernelSerialization::deserialize(kernel_group_opt.value());

            auto aggregator_group_opt = node.getGroup("aggregator");
            if (!aggregator_group_opt) JGAP_LOG_AND_THROW("Missing 'aggregator' group in ManyBodyGapComponent serialization");
            auto aggregator = SerializationRegistry<TransformationAggregator<Dim>>::deserialize(
                aggregator_group_opt.value());

            auto sparse_points = node.loadSparsePoints<Dim>();
            auto coefficients = node.readOptionalDataSet<std::vector<Real>>("coefficients").value_or(std::vector<Real>{});

            return ValuePtr<GapComponent>(ComponentT(aggregator, kernel, sparse_points, coefficients));
        }
    };
}

#endif
