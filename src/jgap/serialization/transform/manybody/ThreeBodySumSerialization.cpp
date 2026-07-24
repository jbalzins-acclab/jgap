#include "ThreeBodySumSerialization.hpp"
#include "jgap/core/atomic/species/composition/Species3AtomicSorted.hpp"
#include "jgap/experimental/transform/manybody/ThreeBodySum.hpp"
#include "jgap/io/log/CurrentLogger.hpp"
#include "jgap/serialization/SerializationNode.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"

namespace jgap {
    template<size_t Dim>
    bool ThreeBodySumSerialization<Dim>::serialize(const ValuePtr<NBodyAggregator<Dim>>& obj,
                                                   SerializationNode& node) const {
        if (auto derived = obj.template as<ThreeBodySum<Dim>>()) {
            node.writeAttribute("name", "ThreeBodySum");
            node.writeAttribute("dim", Dim);
            node.writeAttribute("central_atom_species", derived->getCentralSpecies().symbol());

            auto transformations_group = node.createGroup("transformations");
            int i = 0;
            for (const auto& [species_set, transformation]: derived->getTransformations()) {
                auto transformation_group = transformations_group.createGroup(std::to_string(i++));

                transformation_group.writeAttribute("species_set", species_set.toString());

                auto ct_group = transformation_group.createGroup("cluster_transformation");
                SerializationRegistry<ThreeBodyTransformation<Dim>>::serialize(transformation, ct_group);
            }
            return true;
        }
        return false;
    }

    template<size_t Dim>
    ValuePtr<NBodyAggregator<Dim>> ThreeBodySumSerialization<Dim>::deserialize(const SerializationNode& node) const {
        if (node.readOptionalStringAttribute("name") != "ThreeBodySum" ||
            node.readOptionalSizeAttribute("dim") != Dim) {
            return nullptr;
        }

        auto central_species_symbol = node.readStringAttribute("central_atom_species");
        auto aggregator = std::make_unique<ThreeBodySum<Dim>>(Species(central_species_symbol));

        auto transformations_group_opt = node.getGroup("transformations");
        if (!transformations_group_opt)
            JGAP_LOG_AND_THROW("Missing 'transformations' group in ThreeBodySum serialization");
        const auto& transformations_group = transformations_group_opt.value();

        for (const auto& group_name: transformations_group.getChildNames()) {
            auto transformation_group_opt = transformations_group.getGroup(group_name);
            if (!transformation_group_opt)
                JGAP_LOG_AND_THROW("Missing transformation group in ThreeBodySum serialization");
            const auto& transformation_group = transformation_group_opt.value();

            auto species_encoded = transformation_group.readStringAttribute("species_set");
            Species3AtomicSorted species_set(species_encoded);

            auto ct_group_opt = transformation_group.getGroup("cluster_transformation");
            if (!ct_group_opt)
                JGAP_LOG_AND_THROW("Missing 'cluster_transformation' group in ThreeBodySum serialization");

            auto transformation =
                SerializationRegistry<ThreeBodyTransformation<Dim>>::deserialize(ct_group_opt.value());
            aggregator->extend(species_set, transformation);
        }

        return ValuePtr<NBodyAggregator<Dim>>(std::move(aggregator));
    }

    template class ThreeBodySumSerialization<1>;
    template class ThreeBodySumSerialization<2>;
    template class ThreeBodySumSerialization<3>;
    template class ThreeBodySumSerialization<4>;

    REGISTER_SERIALIZATION(ThreeBodySumSerialization<1>, NBodyAggregator<1>);
    REGISTER_SERIALIZATION(ThreeBodySumSerialization<2>, NBodyAggregator<2>);
    REGISTER_SERIALIZATION(ThreeBodySumSerialization<3>, NBodyAggregator<3>);
    REGISTER_SERIALIZATION(ThreeBodySumSerialization<4>, NBodyAggregator<4>);
}
