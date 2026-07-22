#include "SplinePairTransformationSerialization.hpp"

#include "core/splines/Grid.hpp"
#include "core/splines/HermiteCubicSpline.hpp"
#include "core/transform/nbody/2b/eam/SplinePairTransformation.hpp"
#include "serialization/SerializationNode.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {

    bool SplinePairTransformationSerialization::serialize(const ValuePtr<TwoBodyTransformation<1>>& obj, SerializationNode& node) const {
        if (auto derived = obj.as<SplinePairTransformation>()) {
            node.writeAttribute("name", "SplinePairTransformation");

            const auto& grid = derived->getSpline().getTable();
            node.writeAttribute("dims", std::vector<size_t>(grid.sizes.begin(), grid.sizes.end()));
            node.writeAttribute("spacing", std::vector<Real>(grid.spacing.begin(), grid.spacing.end()));
            node.writeAttribute("origin", std::vector<Real>(grid.origin.begin(), grid.origin.end()));
            node.writeDataSet("data_flat", grid.data_flat);
            return true;
        }
        return false;
    }

    ValuePtr<TwoBodyTransformation<1>> SplinePairTransformationSerialization::deserialize(const SerializationNode& node) const {
        if (node.readOptionalAttribute<std::string>("name") != "SplinePairTransformation") {
            return nullptr;
        }

        auto dims_vec = node.readAttribute<std::vector<size_t>>("dims");
        auto spacing_vec = node.readAttribute<std::vector<Real>>("spacing");
        auto origin_vec = node.readAttribute<std::vector<Real>>("origin");
        auto data_flat = node.readDataSet<std::vector<Real>>("data_flat");

        std::array<size_t, 1> dims = {dims_vec[0]};
        std::array<Real, 1> spacing = {spacing_vec[0]};
        std::array<Real, 1> origin = {origin_vec[0]};

        Grid<1> grid(dims, spacing, origin, data_flat);
        HermiteCubicSpline spline(grid);

        return ValuePtr<TwoBodyTransformation<1>>(SplinePairTransformation(std::move(spline)));
    }

    REGISTER_SERIALIZATION(SplinePairTransformationSerialization, TwoBodyTransformation<1>);
}