#include "TwoBodySumSerialization.hpp"
#include "jgap/core/atomic/species/composition/Species2Atomic.hpp"
#include "jgap/core/transform/manybody/TwoBodySum.hpp"
#include "jgap/io/log/CurrentLogger.hpp"
#include "jgap/serialization/SerializationNode.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"

namespace jgap {
    template<size_t Dim>
    bool TwoBodySumSerialization<Dim>::serialize(const ValuePtr<NBodyAggregator<Dim>>& obj,
                                                 SerializationNode& node) const {
        if (auto derived = obj.template as<TwoBodySum<Dim>>()) {
            node.writeAttribute("name", "TwoBodySum");
            node.writeAttribute("dim", Dim);
            node.writeAttribute("central_atom_species", derived->getCentralSpecies().symbol());

            auto transformations_group = node.createGroup("transformations");
            int i = 0;
            for (const auto& [species_set, transformation]: derived->getTransformations()) {
                auto transformation_group = transformations_group.createGroup(std::to_string(i++));

                transformation_group.writeAttribute("species_set", species_set.toString());

                auto ct_group = transformation_group.createGroup("cluster_transformation");
                SerializationRegistry<TwoBodyTransformation<Dim>>::serialize(transformation, ct_group);
            }
            return true;
        }
        return false;
    }

    template<size_t Dim>
    ValuePtr<NBodyAggregator<Dim>> TwoBodySumSerialization<Dim>::deserialize(const SerializationNode& node) const {
        if (node.readOptionalStringAttribute("name") != "TwoBodySum" ||
            node.readOptionalSizeAttribute("dim") != Dim) {
            return nullptr;
        }

        auto central_species_symbol = node.readStringAttribute("central_atom_species");
        auto aggregator = std::make_unique<TwoBodySum<Dim>>(Species(central_species_symbol));

        auto transformations_group_opt = node.getGroup("transformations");
        if (!transformations_group_opt)
            JGAP_LOG_AND_THROW("Missing 'transformations' group in TwoBodySum serialization");
        const auto& transformations_group = transformations_group_opt.value();

        for (const auto& group_name: transformations_group.getChildNames()) {
            auto transformation_group_opt = transformations_group.getGroup(group_name);
            if (!transformation_group_opt)
                JGAP_LOG_AND_THROW("Missing transformation group in TwoBodySum serialization");
            const auto& transformation_group = transformation_group_opt.value();

            auto species_encoded = transformation_group.readStringAttribute("species_set");
            Species2Atomic species_set(species_encoded);

            auto ct_group_opt = transformation_group.getGroup("cluster_transformation");
            if (!ct_group_opt) JGAP_LOG_AND_THROW("Missing 'cluster_transformation' group in TwoBodySum serialization");

            auto transformation = SerializationRegistry<TwoBodyTransformation<Dim>>::deserialize(ct_group_opt.value());
            aggregator->extend(species_set, transformation);
        }

        return ValuePtr<NBodyAggregator<Dim>>(std::move(aggregator));
    }

    template class TwoBodySumSerialization<1>;
    template class TwoBodySumSerialization<2>;
    template class TwoBodySumSerialization<3>;
    template class TwoBodySumSerialization<4>;

    REGISTER_SERIALIZATION(TwoBodySumSerialization<1>, NBodyAggregator<1>);
    REGISTER_SERIALIZATION(TwoBodySumSerialization<2>, NBodyAggregator<2>);
    REGISTER_SERIALIZATION(TwoBodySumSerialization<3>, NBodyAggregator<3>);
    REGISTER_SERIALIZATION(TwoBodySumSerialization<4>, NBodyAggregator<4>);
}
