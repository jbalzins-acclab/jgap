#ifndef JGAP_TRANSFORMATIONAGGREGATORIMPLSERIALIZATION_HPP
#define JGAP_TRANSFORMATIONAGGREGATORIMPLSERIALIZATION_HPP

#include "core/transform/aggregated/TransformationAggregator.hpp"
#include "core/transform/aggregated/TransformationAggregatorImpl.hpp"
#include "serialization/Serialization.hpp"
#include "serialization/SerializationRegistry.hpp"
#include "core/ValuePtr.hpp"
#include "core/atomic/species/SpeciesSet.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {
    template<size_t Dim, size_t ClusterSize>
    class TransformationAggregatorImplSerialization : public Serialization<TransformationAggregator<Dim>> {
    public:
        bool serialize(const ValuePtr<TransformationAggregator<Dim>>& obj, SerializationNode& node) const override {
            if (auto derived = obj.template as<TransformationAggregatorImpl<Dim, ClusterSize>>()) {
                node.writeAttribute("name", "TransformationAggregatorImpl");
                node.writeAttribute("dim", Dim);
                node.writeAttribute("cluster_size", ClusterSize);
                node.writeAttribute("central_atom_species", derived->getCentralSpecies().symbol());

                auto transformations_group = node.createGroup("transformations");
                int i = 0;
                for (const auto& [species_set, transformation] : derived->getTransformations()) {
                    auto transformation_group = transformations_group.createGroup(std::to_string(i++));

                    std::vector<std::string> species_symbols;
                    species_symbols.push_back(species_set.getRoot().symbol());
                    for (const auto& s : species_set.getNodes()) {
                        species_symbols.push_back(s.symbol());
                    }
                    transformation_group.writeAttribute("species_set", species_symbols);

                    auto ct_group = transformation_group.createGroup("cluster_transformation");
                    SerializationRegistry<ClusterTransformation<Dim, ClusterSize>>::serialize(transformation, ct_group);
                }
                return true;
            }
            return false;
        }

        ValuePtr<TransformationAggregator<Dim>> deserialize(const SerializationNode& node) const override {
            if (node.readOptionalAttribute<std::string>("name") != "TransformationAggregatorImpl" ||
                node.readOptionalAttribute<size_t>("dim") != Dim ||
                node.readOptionalAttribute<size_t>("cluster_size") != ClusterSize) {
                return nullptr;
            }

            auto central_species_symbol = node.readAttribute<std::string>("central_atom_species");
            auto aggregator = std::make_unique<TransformationAggregatorImpl<Dim, ClusterSize>>(Species(central_species_symbol));

            auto transformations_group_opt = node.getGroup("transformations");
            if (!transformations_group_opt) JGAP_LOG_AND_THROW("Missing 'transformations' group in TransformationAggregatorImpl serialization");
            const auto& transformations_group = transformations_group_opt.value();

            for (const auto& group_name : transformations_group.getChildNames()) {
                auto transformation_group_opt = transformations_group.getGroup(group_name);
                if (!transformation_group_opt) JGAP_LOG_AND_THROW("Missing transformation group in TransformationAggregatorImpl serialization");
                const auto& transformation_group = transformation_group_opt.value();

                auto species_symbols = transformation_group.readAttribute<std::vector<std::string>>("species_set");

                Species root(species_symbols[0]);

                if constexpr (ClusterSize == 2) {
                    SpeciesSet<2, HasCentralAtom> species_set(root, Species(species_symbols[1]));
                    auto ct_group_opt = transformation_group.getGroup("cluster_transformation");
                    if (!ct_group_opt) JGAP_LOG_AND_THROW("Missing 'cluster_transformation' group in TransformationAggregatorImpl serialization");
                    auto transformation = SerializationRegistry<ClusterTransformation<Dim, 2>>::deserialize(ct_group_opt.value());
                    aggregator->extend(species_set, transformation);
                } else if constexpr (ClusterSize == 3) {
                    SpeciesSet<3, HasCentralAtom> species_set(root, Species(species_symbols[1]), Species(species_symbols[2]));
                    auto ct_group_opt = transformation_group.getGroup("cluster_transformation");
                    if (!ct_group_opt) JGAP_LOG_AND_THROW("Missing 'cluster_transformation' group in TransformationAggregatorImpl serialization");
                    auto transformation = SerializationRegistry<ClusterTransformation<Dim, 3>>::deserialize(ct_group_opt.value());
                    aggregator->extend(species_set, transformation);
                }
            }

            return ValuePtr<TransformationAggregator<Dim>>(std::move(aggregator));
        }
    };

    // Explicit instantiations
    extern template class TransformationAggregatorImplSerialization<1, 2>;
    extern template class TransformationAggregatorImplSerialization<1, 3>;
    extern template class TransformationAggregatorImplSerialization<2, 2>;
    extern template class TransformationAggregatorImplSerialization<2, 3>;
    extern template class TransformationAggregatorImplSerialization<3, 2>;
    extern template class TransformationAggregatorImplSerialization<3, 3>;
    extern template class TransformationAggregatorImplSerialization<4, 2>;
    extern template class TransformationAggregatorImplSerialization<4, 3>;
}

#endif