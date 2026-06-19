#include <gtest/gtest.h>

#include "core/cutoff/CosCutoff.hpp"
#include "core/splines/Grid.hpp"
#include "core/splines/HermiteCubicSpline.hpp"
#include "core/transform/2b/TwoBodyTransformation.hpp"
#include "core/transform/3b/Angle3bTransformation.hpp"
#include "core/transform/eam/FSGenPairFunction.hpp"
#include "core/transform/eam/CoscutoffPairFunction.hpp"
#include "core/transform/eam/PolycutoffPairFunction.hpp"
#include "core/transform/eam/SplinePairTransformation.hpp"
#include "core/transform/aggregated/TransformationAggregatorImpl.hpp"
#include "serialization/SerializationRegistry.hpp"

using namespace jgap;

namespace {
    std::string tmpFile(const std::string& name) {
        return testing::TempDir() + "/" + name;
    }

    template<typename Base, typename Obj>
    ValuePtr<Base> roundTrip(const Obj& obj, const std::string& file) {
        SerializationRegistry<Base>::serialize(ValuePtr<Base>(obj), file);
        return SerializationRegistry<Base>::deserialize(file);
    }
}

TEST(TestTransformSerialization, TwoBodyTransformation) {
    auto restored = roundTrip<ClusterTransformation<2, 2>>(
        TwoBodyTransformation(CosCutoff(5.0, 1.0)), tmpFile("twobody.h5"));

    auto casted = restored.as<TwoBodyTransformation>();
    ASSERT_NE(casted, nullptr);
    auto cutoff_ptr = casted->getCutoff();
    auto cutoff = cutoff_ptr.as<CosCutoff>();
    ASSERT_NE(cutoff, nullptr);
    EXPECT_DOUBLE_EQ(cutoff->getCutoff(), 5.0);
    EXPECT_DOUBLE_EQ(cutoff->getCutoffTransitionWidth(), 1.0);
}

TEST(TestTransformSerialization, Angle3bTransformation) {
    auto restored = roundTrip<ClusterTransformation<4, 3>>(
        Angle3bTransformation(CosCutoff(4.0, 0.5)), tmpFile("angle3b.h5"));

    auto casted = restored.as<Angle3bTransformation>();
    ASSERT_NE(casted, nullptr);
    auto cutoff_ptr = casted->getCutoff();
    auto cutoff = cutoff_ptr.as<CosCutoff>();
    ASSERT_NE(cutoff, nullptr);
    EXPECT_DOUBLE_EQ(cutoff->getCutoff(), 4.0);
    EXPECT_DOUBLE_EQ(cutoff->getCutoffTransitionWidth(), 0.5);
}

TEST(TestTransformSerialization, FSGenPairFunction) {
    auto restored = roundTrip<ClusterTransformation<1, 2>>(
        FSGenPairFunction(4.5, 3.0, 1.2), tmpFile("fsgen.h5"));

    auto casted = restored.as<FSGenPairFunction>();
    ASSERT_NE(casted, nullptr);
    EXPECT_DOUBLE_EQ(casted->getCutoff(), 4.5);
    EXPECT_DOUBLE_EQ(casted->getDegree(), 3.0);
    EXPECT_DOUBLE_EQ(casted->getPrefactor(), 1.2);
}

TEST(TestTransformSerialization, CoscutoffPairFunction) {
    auto restored = roundTrip<ClusterTransformation<1, 2>>(
        CoscutoffPairFunction(4.0, 1.0, 2.0), tmpFile("coscutoff_pf.h5"));

    auto casted = restored.as<CoscutoffPairFunction>();
    ASSERT_NE(casted, nullptr);
    EXPECT_DOUBLE_EQ(casted->getCutoff(), 4.0);
    EXPECT_DOUBLE_EQ(casted->getRMin(), 1.0);
    EXPECT_DOUBLE_EQ(casted->getPrefactor(), 2.0);
}

TEST(TestTransformSerialization, PolycutoffPairFunction) {
    auto restored = roundTrip<ClusterTransformation<1, 2>>(
        PolycutoffPairFunction(4.0, 1.0, 2.0), tmpFile("polycutoff_pf.h5"));

    auto casted = restored.as<PolycutoffPairFunction>();
    ASSERT_NE(casted, nullptr);
    EXPECT_DOUBLE_EQ(casted->getCutoff(), 4.0);
    EXPECT_DOUBLE_EQ(casted->getRMin(), 1.0);
    EXPECT_DOUBLE_EQ(casted->getPrefactor(), 2.0);
}

TEST(TestTransformSerialization, SplinePairTransformation) {
    Grid<1> grid({4}, {0.5}, {1.0}, {0.0, 1.0, 4.0, 9.0});
    SplinePairTransformation original(HermiteCubicSpline{grid});

    auto restored = roundTrip<ClusterTransformation<1, 2>>(original, tmpFile("spline_pf.h5"));

    auto casted = restored.as<SplinePairTransformation>();
    ASSERT_NE(casted, nullptr);
    const auto& table = casted->getSpline().getTable();
    EXPECT_EQ(table.dims[0], grid.dims[0]);
    EXPECT_DOUBLE_EQ(table.spacing[0], grid.spacing[0]);
    EXPECT_DOUBLE_EQ(table.origin[0], grid.origin[0]);
    ASSERT_EQ(table.data_flat.size(), grid.data_flat.size());
    for (size_t i = 0; i < grid.data_flat.size(); ++i) {
        EXPECT_DOUBLE_EQ(table.data_flat[i], grid.data_flat[i]);
    }
}

TEST(TestTransformSerialization, TransformationAggregatorImpl) {
    TransformationAggregatorImpl<1, 2> aggregator("Fe");
    aggregator.extend({"Fe", "Ni"}, PolycutoffPairFunction(4.0, 1.0, 2.0));

    auto restored = roundTrip<TransformationAggregator<1>>(aggregator, tmpFile("aggregator.h5"));

    auto casted = restored.as<TransformationAggregatorImpl<1, 2>>();
    ASSERT_NE(casted, nullptr);
    EXPECT_EQ(casted->getCentralSpecies().symbol(), "Fe");
    ASSERT_EQ(casted->getTransformations().size(), 1u);

    const auto& [species_set, transformation] = *casted->getTransformations().begin();
    auto pf = transformation.as<PolycutoffPairFunction>();
    ASSERT_NE(pf, nullptr);
    EXPECT_DOUBLE_EQ(pf->getCutoff(), 4.0);
}
