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
        void set_values(Real r, Real val, Real deriv) {
            expected_r = r;
            return_val = val;
            return_deriv = deriv;
        }

        Real evaluate(Real r) const override {
            EXPECT_NEAR(r, expected_r, 1e-9);
            return return_val;
        }

        std::tuple<Real, Real> evaluateAndDifferentiate(Real r) const override {
            EXPECT_NEAR(r, expected_r, 1e-9);
            return {return_val, return_deriv};
        }

        Real getCutoff() const override { return 10.0; }

        MockCutoff* clone() const override { return new MockCutoff(*this); }

    private:
        Real expected_r = 0.0;
        Real return_val = 0.0;
        Real return_deriv = 0.0;
    };
}

TEST(TestPairDistanceTransformation, CorrectlyUsesCutoff) {
    // 2. Setup
    MockCutoff mock_cutoff{};
    Real test_dist = 2.5;
    Real test_val = 0.5;
    Real test_deriv = -0.25;
    mock_cutoff.set_values(test_dist, test_val, test_deriv);

    PairDistanceTransformation trans(mock_cutoff);

    Cluster2 pair;
    pair.separation01.magnitude = test_dist;

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

    EXPECT_EQ(desc_and_derivs.derivatives.size(), 2);
    EXPECT_NEAR(desc_and_derivs.derivatives[0], 1.0, 1e-9); // Derivative of r wrt r is 1
    EXPECT_NEAR(desc_and_derivs.derivatives[1], test_deriv, 1e-9); // Derivative of f_cut wrt r
}
