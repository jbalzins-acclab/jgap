#include <gtest/gtest.h>

#include "core/cutoff/CosCutoff.hpp"
#include "core/kernels/SquaredExpKernel.hpp"
#include "core/potentials/gap/component/GapComponent.hpp"
#include "core/potentials/gap/component/ManyBodyGapComponent.hpp"
#include "core/potentials/gap/component/NBodyGapComponent.hpp"
#include "core/transform/2b/TwoBodyTransformation.hpp"
#include "core/transform/aggregated/TransformationAggregatorImpl.hpp"
#include "core/transform/eam/PolycutoffPairFunction.hpp"
#include "serialization/SerializationRegistry.hpp"

using namespace jgap;

namespace {
    std::string tmpFile(const std::string& name) {
        return testing::TempDir() + "/" + name;
    }

    template<typename Obj>
    ValuePtr<GapComponent> roundTrip(const Obj& component, const std::string& file) {
        SerializationRegistry<GapComponent>::serialize(ValuePtr<GapComponent>(component), file);
        return SerializationRegistry<GapComponent>::deserialize(file);
    }

    void expectCoeffsEq(const std::vector<Real>& a, const std::vector<Real>& b) {
        ASSERT_EQ(a.size(), b.size());
        for (size_t i = 0; i < a.size(); ++i) {
            EXPECT_DOUBLE_EQ(a[i], b[i]);
        }
    }
}

TEST(TestComponentSerialization, NBodyTwoBody) {
    std::vector<Descriptor<2>> sparse_points = {Descriptor<2>{{2.0, 1.0}}, Descriptor<2>{{2.5, 0.5}}};
    std::vector<Real> coefficients = {0.5, -0.3};

    NBodyGapComponent<2, 2, FullSymmetry, SquaredExpKernel<1, 1>> component(
        SpeciesSet<2, FullSymmetry>{"Fe", "Fe"},
        TwoBodyTransformation(CosCutoff(5.0, 1.0)),
        SquaredExpKernel<1, 1>(10.0, {1.3}),
        sparse_points,
        coefficients);

    auto rt = roundTrip(component, tmpFile("nbody_2b.h5"));
    auto restored = rt.as<NBodyGapComponent<2, 2, FullSymmetry, SquaredExpKernel<1, 1>>>();
    ASSERT_NE(restored, nullptr);

    EXPECT_DOUBLE_EQ(restored->getKernel().getEnergyScale(), 10.0);
    EXPECT_DOUBLE_EQ(restored->getKernel().getLengthScales()[0], 1.3);
    expectCoeffsEq(restored->getCoefficients(), coefficients);

    ASSERT_EQ(restored->getSparsePoints().size(), sparse_points.size());
    for (size_t i = 0; i < sparse_points.size(); ++i) {
        EXPECT_DOUBLE_EQ(restored->getSparsePoints()[i].value[0], sparse_points[i].value[0]);
        EXPECT_DOUBLE_EQ(restored->getSparsePoints()[i].value[1], sparse_points[i].value[1]);
    }
}

TEST(TestComponentSerialization, ManyBodyEam) {
    TransformationAggregatorImpl<1, 2> aggregator("Fe");
    aggregator.extend({"Fe", "Ni"}, PolycutoffPairFunction(4.0, 1.0, 2.0));

    std::vector<Descriptor<1>> sparse_points = {Descriptor<1>{{1.5}}};
    std::vector<Real> coefficients = {0.42};

    ManyBodyGapComponent<1, SquaredExpKernel<1, 0>> component(
        aggregator,
        SquaredExpKernel<1, 0>(1.0, {1.0}),
        sparse_points,
        coefficients);

    auto rt = roundTrip(component, tmpFile("manybody_eam.h5"));
    auto restored = rt.as<ManyBodyGapComponent<1, SquaredExpKernel<1, 0>>>();
    ASSERT_NE(restored, nullptr);

    EXPECT_DOUBLE_EQ(restored->getKernel().getEnergyScale(), 1.0);
    expectCoeffsEq(restored->getCoefficients(), coefficients);
    ASSERT_EQ(restored->getSparsePoints().size(), sparse_points.size());
    EXPECT_DOUBLE_EQ(restored->getSparsePoints()[0].value[0], 1.5);

    auto restored_aggregator = restored->getAggregator().as<TransformationAggregatorImpl<1, 2>>();
    ASSERT_NE(restored_aggregator, nullptr);
    EXPECT_EQ(restored_aggregator->getCentralSpecies().symbol(), "Fe");
}
