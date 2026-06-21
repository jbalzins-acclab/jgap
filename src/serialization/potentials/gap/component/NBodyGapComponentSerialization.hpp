#ifndef JGAP_NBODYGAPCOMPONENTSERIALIZATION_HPP
#define JGAP_NBODYGAPCOMPONENTSERIALIZATION_HPP

#include "core/potentials/gap/component/GapComponent.hpp"
#include "core/potentials/gap/component/NBodyGapComponent.hpp"
#include "core/kernels/SquaredExpKernel.hpp"
#include "serialization/Serialization.hpp"
#include "serialization/SerializationRegistry.hpp"
#include "serialization/SerializationNode.hpp"
#include "serialization/kernels/SquaredExpKernelSerialization.hpp"
#include "core/ValuePtr.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {

    // Serializer for the NBodyGapComponent<Dim, ClusterSize, ClusterSym, TKernel> family.
    // The kernel is assumed to be a SquaredExpKernel (the only kernel currently used); TKernel exposes
    // ExpDim/CutoffDim so the matching SquaredExpKernelSerialization can be selected. Explicit
    // instantiations and registrations live in the .cpp.
    template<size_t Dim, size_t ClusterSize, ClusterSymmetry ClusterSym, typename TKernel>
    class NBodyGapComponentSerialization : public Serialization<GapComponent> {
        using ComponentT = NBodyGapComponent<Dim, ClusterSize, ClusterSym, TKernel>;
        using KernelSerialization = SquaredExpKernelSerialization<TKernel::ExpDim, TKernel::CutoffDim>;

    public:
        bool serialize(const ValuePtr<GapComponent>& obj, SerializationNode& node) const override {
            auto derived = obj.template as<ComponentT>();
            if (!derived) {
                return false;
            }

            node.writeAttribute("name", "NBodyGapComponent");
            node.writeAttribute("dim", Dim);
            node.writeAttribute("cluster_size", ClusterSize);
            node.writeAttribute("cluster_symmetry", symmetryName());

            // Species set: for Symmetric this is just the nodes, for HasCentralAtom the root comes first.
            std::vector<std::string> species_symbols;
            if constexpr (ClusterSym == HasCentralAtom) {
                species_symbols.push_back(derived->getSpecies().getRoot().symbol());
            }
            for (const auto& s : derived->getSpecies().getNodes()) {
                species_symbols.push_back(s.symbol());
            }
            node.writeAttribute("species_set", species_symbols);

            auto kernel_group = node.createGroup("kernel");
            KernelSerialization::serialize(derived->getKernel(), kernel_group);

            auto transformation_group = node.createGroup("transformation");
            SerializationRegistry<ClusterTransformation<Dim, ClusterSize>>::serialize(
                derived->getTransformation(), transformation_group);

            node.saveSparsePoints<Dim>(derived->getSparsePoints());

            const auto& coefficients = derived->getCoefficients();
            if (!coefficients.empty()) {
                node.writeDataSet("coefficients", coefficients);
            }

            return true;
        }

        ValuePtr<GapComponent> deserialize(const SerializationNode& node) const override {
            if (node.readOptionalAttribute<std::string>("name") != "NBodyGapComponent" ||
                node.readOptionalAttribute<size_t>("dim") != Dim ||
                node.readOptionalAttribute<size_t>("cluster_size") != ClusterSize ||
                node.readOptionalAttribute<std::string>("cluster_symmetry") != symmetryName()) {
                return nullptr;
            }

            auto species_symbols = node.readAttribute<std::vector<std::string>>("species_set");

            auto kernel_group_opt = node.getGroup("kernel");
            if (!kernel_group_opt) JGAP_LOG_AND_THROW("Missing 'kernel' group in NBodyGapComponent serialization");
            TKernel kernel = KernelSerialization::deserialize(kernel_group_opt.value());

            auto transformation_group_opt = node.getGroup("transformation");
            if (!transformation_group_opt) JGAP_LOG_AND_THROW("Missing 'transformation' group in NBodyGapComponent serialization");
            auto transformation = SerializationRegistry<ClusterTransformation<Dim, ClusterSize>>::deserialize(
                transformation_group_opt.value());

            auto sparse_points = node.loadSparsePoints<Dim>();
            auto coefficients = node.readOptionalDataSet<std::vector<Real>>("coefficients").value_or(std::vector<Real>{});

            return ValuePtr<GapComponent>(
                ComponentT(makeSpeciesSet(species_symbols), transformation, kernel, sparse_points, coefficients));
        }

    private:
        static std::string symmetryName() {
            return ClusterSym == HasCentralAtom ? "HasCentralAtom" : "Symmetric";
        }

        static SpeciesSet<ClusterSize, ClusterSym> makeSpeciesSet(const std::vector<std::string>& symbols) {
            if constexpr (ClusterSym == FullSymmetry && ClusterSize == 2) {
                if (symbols.size() != 2) JGAP_LOG_AND_THROW("Expected 2 species for NBodyGapComponent<*,2,Symmetric>");
                return SpeciesSet<2, FullSymmetry>(Species(symbols[0]), Species(symbols[1]));
            } else if constexpr (ClusterSym == HasCentralAtom && ClusterSize == 3) {
                if (symbols.size() != 3) JGAP_LOG_AND_THROW("Expected 3 species for NBodyGapComponent<*,3,HasCentralAtom>");
                return SpeciesSet<3, HasCentralAtom>(Species(symbols[0]), Species(symbols[1]), Species(symbols[2]));
            } else {
                JGAP_LOG_AND_THROW("Species set reconstruction not implemented for this NBodyGapComponent specialization");
            }
        }
    };
}

#endif
