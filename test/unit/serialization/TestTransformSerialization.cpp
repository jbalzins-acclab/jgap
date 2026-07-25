#include <gtest/gtest.h>

#include "jgap/core/cutoff/CosCutoff.hpp"
#include "jgap/core/splines/Grid.hpp"
#include "jgap/core/splines/HermiteCubicSpline.hpp"
#include "jgap/core/transform/manybody/TwoBodySum.hpp"
#include "jgap/core/transform/nbody/2b/PairDistanceTransformation.hpp"
#include "jgap/core/transform/nbody/2b/eam/CoscutoffPairFunction.hpp"
#include "jgap/core/transform/nbody/2b/eam/FSGenPairFunction.hpp"
#include "jgap/core/transform/nbody/2b/eam/PolycutoffPairFunction.hpp"
#include "jgap/core/transform/nbody/2b/eam/SplinePairTransformation.hpp"
#include "jgap/core/transform/nbody/3b/Angle3bTransformation.hpp"
#include "jgap/experimental/transform/nbody/3b/CutoffJK3bTransformation.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"

using namespace jgap;

namespace {
    std::string tmpFile(const std::string& name) { return testing::TempDir() + "/" + name; }

    template<typename Base, typename Obj>
    ValuePtr<Base> roundTrip(const Obj& obj, const std::string& file) {
        SerializationRegistry<Base>::serialize(ValuePtr<Base>(obj), file);
        return SerializationRegistry<Base>::deserialize(file);
    }
}

TEST(TestTransformSerialization, TwoPairDistanceTransformation) {
    auto restored =
        roundTrip<TwoBodyTransformation<2>>(PairDistanceTransformation(CosCutoff(5.0, 1.0)), tmpFile("twobody.h5"));

    auto casted = restored.as<PairDistanceTransformation>();
    ASSERT_NE(casted, nullptr);
    auto cutoff_ptr = casted->getCutoffFunction();
    auto cutoff = cutoff_ptr.as<CosCutoff>();
    ASSERT_NE(cutoff, nullptr);
    EXPECT_DOUBLE_EQ(cutoff->getCutoff(), 5.0);
    EXPECT_DOUBLE_EQ(cutoff->getCutoffTransitionWidth(), 1.0);
}

TEST(TestTransformSerialization, Angle3bTransformation) {
    auto restored =
        roundTrip<ThreeBodyTransformation<4>>(Angle3bTransformation(CosCutoff(4.0, 0.5)), tmpFile("angle3b.h5"));

    auto casted = restored.as<Angle3bTransformation>();
    ASSERT_NE(casted, nullptr);
    auto cutoff_ptr = casted->getCutoffFunction();
    auto cutoff = cutoff_ptr.as<CosCutoff>();
    ASSERT_NE(cutoff, nullptr);
    EXPECT_DOUBLE_EQ(cutoff->getCutoff(), 4.0);
    EXPECT_DOUBLE_EQ(cutoff->getCutoffTransitionWidth(), 0.5);
}

TEST(TestTransformSerialization, CutoffJK3bTransformation) {
    auto restored = roundTrip<ThreeBodyTransformation<4>>(
        CutoffJK3bTransformation(CosCutoff(4.0, 0.5), CosCutoff(3.0, 0.2)), tmpFile("cutoff_jk3b.h5")
    );

    auto casted = restored.as<CutoffJK3bTransformation>();
    ASSERT_NE(casted, nullptr);

    auto main_cut = casted->getMainCutoffFunction().as<CosCutoff>();
    ASSERT_NE(main_cut, nullptr);
    EXPECT_DOUBLE_EQ(main_cut->getCutoff(), 4.0);

    auto cut_12 = casted->getCutoff12Function().as<CosCutoff>();
    ASSERT_NE(cut_12, nullptr);
    EXPECT_DOUBLE_EQ(cut_12->getCutoff(), 3.0);
}

TEST(TestTransformSerialization, FSGenPairFunction) {
    auto restored = roundTrip<TwoBodyTransformation<1>>(FSGenPairFunction(4.5, 3.0, 1.2), tmpFile("fsgen.h5"));

    auto casted = restored.as<FSGenPairFunction>();
    ASSERT_NE(casted, nullptr);
    EXPECT_DOUBLE_EQ(casted->getCutoff(), 4.5);
    EXPECT_DOUBLE_EQ(casted->getDegree(), 3.0);
    EXPECT_DOUBLE_EQ(casted->getPrefactor(), 1.2);
}

TEST(TestTransformSerialization, CoscutoffPairFunction) {
    auto restored =
        roundTrip<TwoBodyTransformation<1>>(CoscutoffPairFunction(4.0, 1.0, 2.0), tmpFile("coscutoff_pf.h5"));

    auto casted = restored.as<CoscutoffPairFunction>();
    ASSERT_NE(casted, nullptr);
    EXPECT_DOUBLE_EQ(casted->getCutoff(), 4.0);
    EXPECT_DOUBLE_EQ(casted->getRMin(), 1.0);
    EXPECT_DOUBLE_EQ(casted->getPrefactor(), 2.0);
}

TEST(TestTransformSerialization, PolycutoffPairFunction) {
    auto restored =
        roundTrip<TwoBodyTransformation<1>>(PolycutoffPairFunction(4.0, 1.0, 2.0), tmpFile("polycutoff_pf.h5"));

    auto casted = restored.as<PolycutoffPairFunction>();
    ASSERT_NE(casted, nullptr);
    EXPECT_DOUBLE_EQ(casted->getCutoff(), 4.0);
    EXPECT_DOUBLE_EQ(casted->getRMin(), 1.0);
    EXPECT_DOUBLE_EQ(casted->getPrefactor(), 2.0);
}

TEST(TestTransformSerialization, SplinePairTransformation) {
    Grid<1> grid({4}, {0.5}, {1.0}, {0.0, 1.0, 4.0, 9.0});
    SplinePairTransformation original(HermiteCubicSpline{grid});

    auto restored = roundTrip<TwoBodyTransformation<1>>(original, tmpFile("spline_pf.h5"));

    auto casted = restored.as<SplinePairTransformation>();
    ASSERT_NE(casted, nullptr);
    const auto& table = casted->getSpline().getTable();
    EXPECT_EQ(table.sizes[0], grid.sizes[0]);
    EXPECT_DOUBLE_EQ(table.spacing[0], grid.spacing[0]);
    EXPECT_DOUBLE_EQ(table.origin[0], grid.origin[0]);
    ASSERT_EQ(table.data_flat.size(), grid.data_flat.size());
    for (size_t i = 0; i < grid.data_flat.size(); ++i) {
        EXPECT_DOUBLE_EQ(table.data_flat[i], grid.data_flat[i]);
    }
}

TEST(TestTransformSerialization, NBodyAggregatorImpl) {
    TwoBodySum<1> aggregator("Fe");
    aggregator.extend({"Fe", "Ni"}, PolycutoffPairFunction(4.0, 1.0, 2.0));

    auto restored = roundTrip<NBodyAggregator<1>>(aggregator, tmpFile("aggregator.h5"));

    auto casted = restored.as<TwoBodySum<1>>();
    ASSERT_NE(casted, nullptr);
    EXPECT_EQ(casted->getCentralSpecies().symbol(), "Fe");
    ASSERT_EQ(casted->getTransformations().size(), 1u);

    const auto& [species_set, transformation] = *casted->getTransformations().begin();
    auto pf = transformation.as<PolycutoffPairFunction>();
    ASSERT_NE(pf, nullptr);
    EXPECT_DOUBLE_EQ(pf->getCutoff(), 4.0);
}
