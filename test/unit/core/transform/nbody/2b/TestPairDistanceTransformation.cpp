#include <gtest/gtest.h>
#include <memory>
#include "jgap/core/atomic/geometry/Cluster2.hpp"
#include "jgap/core/cutoff/CutoffFunction.hpp"
#include "jgap/core/transform/nbody/2b/PairDistanceTransformation.hpp"

using namespace jgap;

namespace {
    // 1. Mock Cutoff Function for complete control
    class MockCutoff : public CutoffFunction {
    public:
        // Set the values that the mock will return
        void set_values(double r, double val, double deriv) {
            expected_r = r;
            return_val = val;
            return_deriv = deriv;
        }

        double evaluate(double r) const override {
            EXPECT_NEAR(r, expected_r, 1e-9);
            return return_val;
        }

        std::tuple<double, double> evaluateAndDifferentiate(double r) const override {
            EXPECT_NEAR(r, expected_r, 1e-9);
            return {return_val, return_deriv};
        }

        double getCutoff() const override { return 10.0; }

        MockCutoff* clone() const override { return new MockCutoff(*this); }

    private:
        double expected_r = 0.0;
        double return_val = 0.0;
        double return_deriv = 0.0;
    };
}

TEST(TestPairDistanceTransformation, CorrectlyUsesCutoff) {
    // 2. Setup
    MockCutoff mock_cutoff{};
    double test_dist = 2.5;
    double test_val = 0.5;
    double test_deriv = -0.25;
    mock_cutoff.set_values(test_dist, test_val, test_deriv);

    PairDistanceTransformation trans(mock_cutoff);

    Cluster2 pair;
    pair.separation01.magnitude = test_dist;
    pair.separation01.direction = Vector3{1.0, 0.0, 0.0};

    // 3. Test evaluate()
    auto desc = trans.evaluate(pair);
    EXPECT_EQ(desc.size(), 2);
    EXPECT_NEAR(desc[0], test_dist, 1e-9); // Should return the distance
    EXPECT_NEAR(desc[1], test_val, 1e-9); // Should return the mocked cutoff value

    // 4. Test evaluateAndDifferentiate()
    auto desc_and_derivs = trans.evaluateAndDifferentiate(pair);
    EXPECT_EQ(desc_and_derivs.value.size(), 2);
    EXPECT_NEAR(desc_and_derivs.value[0], test_dist, 1e-9);
    EXPECT_NEAR(desc_and_derivs.value[1], test_val, 1e-9);

    EXPECT_EQ(desc_and_derivs.grad_r1.size(), 2);
    EXPECT_NEAR(desc_and_derivs.grad_r1[0].x, 1.0, 1e-9);
    EXPECT_NEAR(desc_and_derivs.grad_r1[0].y, 0.0, 1e-9);
    EXPECT_NEAR(desc_and_derivs.grad_r1[0].z, 0.0, 1e-9);

    EXPECT_NEAR(desc_and_derivs.grad_r1[1].x, test_deriv, 1e-9);
    EXPECT_NEAR(desc_and_derivs.grad_r1[1].y, 0.0, 1e-9);
    EXPECT_NEAR(desc_and_derivs.grad_r1[1].z, 0.0, 1e-9);
}
