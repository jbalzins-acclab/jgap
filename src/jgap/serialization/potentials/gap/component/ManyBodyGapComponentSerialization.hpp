#ifndef JGAP_MANYBODYGAPCOMPONENTSERIALIZATION_HPP
#define JGAP_MANYBODYGAPCOMPONENTSERIALIZATION_HPP

#include <vector>
#include "jgap/core/ValuePtr.hpp"
#include "jgap/core/atomic/descriptor/Descriptor.hpp"
#include "jgap/core/kernels/SquaredExpKernel.hpp"
#include "jgap/core/potentials/gap/component/GapComponent.hpp"
#include "jgap/core/potentials/gap/component/ManyBodyGapComponent.hpp"
#include "jgap/io/log/CurrentLogger.hpp"
#include "jgap/serialization/Serialization.hpp"
#include "jgap/serialization/SerializationNode.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"
#include "jgap/core/kernels/Kernel.hpp"
#include "jgap/serialization/Serialization.hpp"
#include "jgap/serialization/SerializationNode.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"

namespace jgap {

    // Serializer for the ManyBodyGapComponent<Dim, TKernel> family.
    template<size_t Dim, typename TKernel>
    class ManyBodyGapComponentSerialization : public Serialization<GapComponent> {
        using ComponentT = ManyBodyGapComponent<Dim, TKernel>;

    public:
        bool serialize(const ValuePtr<GapComponent>& obj, SerializationNode& node) const override {
            const auto derived = obj.as<ComponentT>();
            if (!derived) {
                return false;
            }

            node.writeAttribute("name", "ManyBodyGapComponent");
            node.writeAttribute("dim", Dim);

            auto kernel_group = node.createGroup("kernel");
            ValuePtr<Kernel<Dim>> kernel_ptr = derived->getKernel();
            SerializationRegistry<Kernel<Dim>>::serialize(kernel_ptr, kernel_group);

            auto aggregator_group = node.createGroup("aggregator");
            SerializationRegistry<NBodyAggregator<Dim>>::serialize(derived->getAggregator(), aggregator_group);

            node.writeDataSet("sparse_points", derived->getSparsePoints());

            const auto& coefficients = derived->getCoefficients();
            if (!coefficients.empty()) {
                node.writeDataSet("coefficients", coefficients);
            }

            return true;
        }

        ValuePtr<GapComponent> deserialize(const SerializationNode& node) const override {
            if (node.readOptionalStringAttribute("name") != "ManyBodyGapComponent" ||
                node.readOptionalSizeAttribute("dim") != Dim) {
                return nullptr;
            }

            auto kernel_group_opt = node.getGroup("kernel");
            if (!kernel_group_opt) JGAP_LOG_AND_THROW("Missing 'kernel' group in ManyBodyGapComponent serialization");
            auto kernel_ptr =
                SerializationRegistry<Kernel<Dim>>::deserialize(kernel_group_opt.value());
            auto* kernel_typed = kernel_ptr.template as<TKernel>();
            if (!kernel_typed) {
                JGAP_LOG_AND_THROW("Deserialized kernel is not of expected type {}", typeName<TKernel>());
            }
            TKernel kernel = *kernel_typed;

            auto aggregator_group_opt = node.getGroup("aggregator");
            if (!aggregator_group_opt)
                JGAP_LOG_AND_THROW("Missing 'aggregator' group in ManyBodyGapComponent serialization");
            auto aggregator = SerializationRegistry<NBodyAggregator<Dim>>::deserialize(aggregator_group_opt.value());

            auto sparse_points = node.readDescriptors<Dim>("sparse_points");
            auto coefficients =
                node.readOptionalRealVectorDataSet("coefficients").value_or(std::vector<Real>{});

            return ValuePtr<GapComponent>(ComponentT(aggregator, kernel, sparse_points, coefficients));
        }
    };
}

#endif
