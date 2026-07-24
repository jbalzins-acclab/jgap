#include <gtest/gtest.h>

#include "jgap/core/atomic/species/composition/Species3AtomicSorted.hpp"
#include "jgap/core/cutoff/CosCutoff.hpp"
#include "jgap/core/kernels/SquaredExpKernel.hpp"
#include "jgap/core/potentials/gap/component/AtomicThreeBodyGapComponent.hpp"
#include "jgap/core/potentials/gap/component/ManyBodyGapComponent.hpp"
#include "jgap/core/potentials/gap/component/TwoBodyGapComponent.hpp"
#include "jgap/core/transform/manybody/TwoBodySum.hpp"
#include "jgap/core/transform/nbody/2b/PairDistanceTransformation.hpp"
#include "jgap/core/transform/nbody/2b/eam/PolycutoffPairFunction.hpp"
#include "jgap/core/transform/nbody/3b/Angle3bTransformation.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"

using namespace jgap;

namespace {
    std::string tmpFile(const std::string& name) { return testing::TempDir() + "/" + name; }

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

    auto transform_ptr = ValuePtr<TwoBodyTransformation<2>>(PairDistanceTransformation(CosCutoff(5.0, 1.0)));
    TwoBodyGapComponent<2, SquaredExpKernel<1, 1>> component(Species2Sorted(Species("Fe"), Species("Fe")),
                                                             transform_ptr, SquaredExpKernel<1, 1>(10.0, {1.3}),
                                                             sparse_points, coefficients);

    auto rt = roundTrip(component, tmpFile("nbody_2b.h5"));
    auto restored = rt.as<TwoBodyGapComponent<2, SquaredExpKernel<1, 1>>>();
    ASSERT_NE(restored, nullptr);

    EXPECT_DOUBLE_EQ(restored->getKernel().getEnergyScale(), 10.0);
    EXPECT_DOUBLE_EQ(restored->getKernel().getLengthScales()[0], 1.3);
    expectCoeffsEq(restored->getCoefficients(), coefficients);

    ASSERT_EQ(restored->getSparsePoints().size(), sparse_points.size());
    for (size_t i = 0; i < sparse_points.size(); ++i) {
        EXPECT_DOUBLE_EQ(restored->getSparsePoints()[i][0], sparse_points[i][0]);
        EXPECT_DOUBLE_EQ(restored->getSparsePoints()[i][1], sparse_points[i][1]);
    }
}

TEST(TestComponentSerialization, ManyBodyEam) {
    TwoBodySum<1> aggregator("Fe");
    aggregator.extend(Species2Atomic("Fe|Ni"),
                      ValuePtr<TwoBodyTransformation<1>>(PolycutoffPairFunction(4.0, 1.0, 2.0)));

    std::vector<Descriptor<1>> sparse_points = {Descriptor<1>{{1.5}}};
    std::vector<Real> coefficients = {0.42};

    ManyBodyGapComponent<1, SquaredExpKernel<1, 0>> component(aggregator, SquaredExpKernel<1, 0>(1.0, {1.0}),
                                                              sparse_points, coefficients);

    auto rt = roundTrip(component, tmpFile("manybody_eam.h5"));
    auto restored = rt.as<ManyBodyGapComponent<1, SquaredExpKernel<1, 0>>>();
    ASSERT_NE(restored, nullptr);

    EXPECT_DOUBLE_EQ(restored->getKernel().getEnergyScale(), 1.0);
    expectCoeffsEq(restored->getCoefficients(), coefficients);
    ASSERT_EQ(restored->getSparsePoints().size(), sparse_points.size());
    EXPECT_DOUBLE_EQ(restored->getSparsePoints()[0][0], 1.5);

    auto restored_aggregator = restored->getAggregator().as<TwoBodySum<1>>();
    ASSERT_NE(restored_aggregator, nullptr);
    EXPECT_EQ(restored_aggregator->getCentralSpecies().symbol(), "Fe");
}

TEST(TestComponentSerialization, AtomicThreeBody) {
    std::vector<Descriptor<4>> sparse_points = {Descriptor<4>{{6.0, 0.0, 3.0, 1.0}},
                                                Descriptor<4>{{5.0, 0.5, 2.5, 0.8}}};
    std::vector<Real> coefficients = {0.15, -0.25};

    auto transform_ptr = ValuePtr<ThreeBodyTransformation<4>>(Angle3bTransformation(CosCutoff(6.0, 1.0)));
    AtomicThreeBodyGapComponent<4, SquaredExpKernel<3, 1>> component(Species3AtomicSorted("Fe|Fe,Ni"), transform_ptr,
                                                                     SquaredExpKernel<3, 1>(5.0, {1.0, 1.5, 2.0}),
                                                                     sparse_points, coefficients);

    auto rt = roundTrip(component, tmpFile("atomic_3b.h5"));
    auto restored = rt.as<AtomicThreeBodyGapComponent<4, SquaredExpKernel<3, 1>>>();
    ASSERT_NE(restored, nullptr);

    EXPECT_DOUBLE_EQ(restored->getKernel().getEnergyScale(), 5.0);
    EXPECT_DOUBLE_EQ(restored->getKernel().getLengthScales()[0], 1.0);
    EXPECT_DOUBLE_EQ(restored->getKernel().getLengthScales()[1], 1.5);
    EXPECT_DOUBLE_EQ(restored->getKernel().getLengthScales()[2], 2.0);
    expectCoeffsEq(restored->getCoefficients(), coefficients);

    ASSERT_EQ(restored->getSparsePoints().size(), sparse_points.size());
    for (size_t i = 0; i < sparse_points.size(); ++i) {
        for (size_t d = 0; d < 4; ++d) {
            EXPECT_DOUBLE_EQ(restored->getSparsePoints()[i][d], sparse_points[i][d]);
        }
    }
}
