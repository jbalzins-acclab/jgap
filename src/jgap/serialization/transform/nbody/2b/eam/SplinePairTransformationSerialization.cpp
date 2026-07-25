#include "SplinePairTransformationSerialization.hpp"

#include "jgap/core/splines/Grid.hpp"
#include "jgap/core/splines/HermiteCubicSpline.hpp"
#include "jgap/core/transform/nbody/2b/eam/SplinePairTransformation.hpp"
#include "jgap/serialization/SerializationNode.hpp"
#include "jgap/io/log/CurrentLogger.hpp"

namespace jgap {

    bool SplinePairTransformationSerialization::serialize(const ValuePtr<TwoBodyTransformation<1>>& obj, SerializationNode& node) const {
        if (auto derived = obj.as<SplinePairTransformation>()) {
            node.writeAttribute("name", "SplinePairTransformation");

            const auto& grid = derived->getSpline().getTable();
            node.writeDataSet("dims", std::vector<size_t>(grid.sizes.begin(), grid.sizes.end()));
            node.writeDataSet("spacing", std::vector<Real>(grid.spacing.begin(), grid.spacing.end()));
            node.writeDataSet("origin", std::vector<Real>(grid.origin.begin(), grid.origin.end()));
            node.writeDataSet("data_flat", grid.data_flat);
            return true;
        }
        return false;
    }

    ValuePtr<TwoBodyTransformation<1>> SplinePairTransformationSerialization::deserialize(const SerializationNode& node) const {
        if (node.readOptionalStringAttribute("name") != "SplinePairTransformation") {
            return nullptr;
        }

        auto dims_vec = node.readSizeVectorDataSet("dims");
        auto spacing_vec = node.readRealVectorDataSet("spacing");
        auto origin_vec = node.readRealVectorDataSet("origin");
        auto data_flat = node.readRealVectorDataSet("data_flat");

        std::array<size_t, 1> dims = {dims_vec[0]};
        std::array<Real, 1> spacing = {spacing_vec[0]};
        std::array<Real, 1> origin = {origin_vec[0]};

        Grid<1> grid(dims, spacing, origin, data_flat);
        HermiteCubicSpline spline(grid);

        return ValuePtr<TwoBodyTransformation<1>>(SplinePairTransformation(std::move(spline)));
    }

    REGISTER_SERIALIZATION(SplinePairTransformationSerialization, TwoBodyTransformation<1>);
}