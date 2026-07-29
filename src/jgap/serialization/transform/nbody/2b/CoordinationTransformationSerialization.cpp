#include "CoordinationTransformationSerialization.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"
#include <string>

namespace jgap {

    template<size_t Dim>
    bool CoordinationTransformationSerialization<Dim>::serialize(const ValuePtr<TwoBodyTransformation<Dim>>& base_obj, SerializationNode& node) const {
        if (auto obj = dynamic_cast<const CoordinationTransformation<Dim>*>(base_obj.get())) {
            node.writeAttribute("name", "CoordinationTransformation");
            node.writeAttribute("dimensions", Dim);

            const auto& ranges = obj->getRanges();
            std::vector<Real> r_mins;
            std::vector<Real> r_maxs;
            for (size_t i = 0; i < Dim; ++i) {
                r_mins.push_back(ranges[i].first);
                r_maxs.push_back(ranges[i].second);
            }

            node.writeDataSet("r_mins", r_mins);
            node.writeDataSet("r_maxs", r_maxs);

            return true;
        }
        return false;
    }

    template<size_t Dim>
    ValuePtr<TwoBodyTransformation<Dim>> CoordinationTransformationSerialization<Dim>::deserialize(const SerializationNode& node) const {
        if (node.readOptionalStringAttribute("name") != "CoordinationTransformation") {
            return nullptr;
        }

        auto serialized_dim = node.readSizeAttribute("dimensions");
        if (serialized_dim != Dim) {
            JGAP_LOG_AND_THROW("Dimension mismatch: expected {}, got {}", Dim, serialized_dim);
        }

        auto r_mins = node.readRealVectorDataSet("r_mins");
        auto r_maxs = node.readRealVectorDataSet("r_maxs");

        if (r_mins.size() != Dim || r_maxs.size() != Dim) {
            JGAP_LOG_AND_THROW("Invalid array sizes for r_mins or r_maxs in CoordinationTransformation");
        }

        std::array<std::pair<Real, Real>, Dim> ranges;
        for (size_t i = 0; i < Dim; ++i) {
            ranges[i] = {r_mins[i], r_maxs[i]};
        }

        return CoordinationTransformation<Dim>(ranges);
    }

    template class CoordinationTransformationSerialization<1>;
    template class CoordinationTransformationSerialization<2>;
    template class CoordinationTransformationSerialization<3>;
    template class CoordinationTransformationSerialization<4>;

    using CoordinationSerialization1 = CoordinationTransformationSerialization<1>;
    using CoordinationSerialization2 = CoordinationTransformationSerialization<2>;
    using CoordinationSerialization3 = CoordinationTransformationSerialization<3>;
    using CoordinationSerialization4 = CoordinationTransformationSerialization<4>;

    REGISTER_SERIALIZATION(CoordinationSerialization1, TwoBodyTransformation<1>);
    REGISTER_SERIALIZATION(CoordinationSerialization2, TwoBodyTransformation<2>);
    REGISTER_SERIALIZATION(CoordinationSerialization3, TwoBodyTransformation<3>);
    REGISTER_SERIALIZATION(CoordinationSerialization4, TwoBodyTransformation<4>);
}
