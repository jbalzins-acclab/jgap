#ifndef JGAP_TWOBODYGAPCOMPONENTSERIALIZATION_HPP
#define JGAP_TWOBODYGAPCOMPONENTSERIALIZATION_HPP

#include <vector>
#include "core/ValuePtr.hpp"
#include "core/atomic/descriptor/Descriptor.hpp"
#include "core/kernels/SquaredExpKernel.hpp"
#include "core/potentials/gap/component/GapComponent.hpp"
#include "core/potentials/gap/component/TwoBodyGapComponent.hpp"
#include "io/log/CurrentLogger.hpp"
#include "serialization/Serialization.hpp"
#include "serialization/SerializationNode.hpp"
#include "serialization/SerializationRegistry.hpp"
#include "serialization/kernels/SquaredExpKernelSerialization.hpp"

namespace jgap {

    template<size_t Dim, typename TKernel>
    class TwoBodyGapComponentSerialization : public Serialization<GapComponent> {
        using ComponentT = TwoBodyGapComponent<Dim, TKernel>;
        using KernelSerialization = SquaredExpKernelSerialization<TKernel::ExpDim, TKernel::CutoffDim>;

    public:
        bool serialize(const ValuePtr<GapComponent>& obj, SerializationNode& node) const override {
            const auto derived = obj.template as<ComponentT>();
            if (!derived) {
                return false;
            }

            node.writeAttribute("name", "TwoBodyGapComponent");
            node.writeAttribute("dim", Dim);

            node.writeAttribute("species_set", derived->getSpecies().toString());

            auto kernel_group = node.createGroup("kernel");
            KernelSerialization::serialize(derived->getKernel(), kernel_group);

            auto transformation_group = node.createGroup("transformation");
            SerializationRegistry<TwoBodyTransformation<Dim>>::serialize(derived->getTransformation(),
                                                                         transformation_group);

            node.writeDataSet("sparse_points", derived->getSparsePoints());

            const auto& coefficients = derived->getCoefficients();
            if (!coefficients.empty()) {
                node.writeDataSet("coefficients", coefficients);
            }

            return true;
        }

        ValuePtr<GapComponent> deserialize(const SerializationNode& node) const override {
            if (node.readOptionalAttribute<std::string>("name") != "TwoBodyGapComponent" ||
                node.readOptionalAttribute<size_t>("dim") != Dim) {
                return nullptr;
            }

            auto species_encoded = node.readAttribute<std::string>("species_set");
            Species2Sorted species_set(species_encoded);

            auto kernel_group_opt = node.getGroup("kernel");
            if (!kernel_group_opt) JGAP_LOG_AND_THROW("Missing 'kernel' group in TwoBodyGapComponent serialization");
            TKernel kernel = KernelSerialization::deserialize(kernel_group_opt.value());

            auto transformation_group_opt = node.getGroup("transformation");
            if (!transformation_group_opt)
                JGAP_LOG_AND_THROW("Missing 'transformation' group in TwoBodyGapComponent serialization");
            auto transformation =
                SerializationRegistry<TwoBodyTransformation<Dim>>::deserialize(transformation_group_opt.value());

            auto sparse_points = node.readDataSet<std::vector<Descriptor<Dim>>>("sparse_points");
            auto coefficients =
                node.readOptionalDataSet<std::vector<Real>>("coefficients").value_or(std::vector<Real>{});

            return ValuePtr<GapComponent>(ComponentT(species_set, transformation, kernel, sparse_points, coefficients));
        }
    };
}

#endif
