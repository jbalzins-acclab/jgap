#include <gtest/gtest.h>
#include <map>
#include <tuple>
#include "core/cutoff/CutoffFunction.hpp"
#include "core/transform/nbody/3b/Angle3bTransformation.hpp"

using namespace jgap;

namespace {
    // 1. Mock Cutoff Function for complete control
    class MockCutoff : public CutoffFunction {
    public:
        // Set the values that the mock will return for a specific distance
        void add_expected_call(Real r, Real val, Real deriv) { expected_calls[r] = {val, deriv}; }

        Real evaluate(Real r) const override {
            auto it = expected_calls.find(r);
            if (it != expected_calls.end()) {
                return it->second.first;
            }
            ADD_FAILURE() << "evaluate called with unexpected r: " << r;
            return 0.0;
        }

        std::tuple<Real, Real> evaluateAndDifferentiate(Real r) const override {
            auto it = expected_calls.find(r);
            if (it != expected_calls.end()) {
                return it->second;
            }
            ADD_FAILURE() << "evaluateAndDifferentiate called with unexpected r: " << r;
            return {0.0, 0.0};
        }

        Real getCutoff() const override { return 10.0; }

        MockCutoff* clone() const override { return new MockCutoff(*this); }

    private:
        // Using a map is fine as long as the keys are precise and the header is included.
        std::map<Real, std::pair<Real, Real>> expected_calls;
    };
}

TEST(TestAngle3bTransformation, CorrectlyUsesCutoff) {
    // 2. Setup
    MockCutoff mock_cutoff{};

    // Define the properties of our test cluster
    Real r01 = 2.0, r02 = 3.0, r12 = 4.0;

    // Program the mock to return specific values for the distances it will be called with
    Real val01 = 0.8, deriv01 = -0.1;
    Real val02 = 0.6, deriv02 = -0.2;
    mock_cutoff.add_expected_call(r01, val01, deriv01);
    mock_cutoff.add_expected_call(r02, val02, deriv02);

    Angle3bTransformation trans(mock_cutoff);

    Cluster3 triplet;
    triplet.separationBetween(0, 1) = r01;
    triplet.separationBetween(0, 2) = r02;
    triplet.separationBetween(1, 2) = r12;

    // 3. Test evaluate()
    auto desc = trans.evaluate(triplet);
    EXPECT_EQ(desc.size(), 4);
    EXPECT_NEAR(desc[0], r01 + r02, 1e-9);
    EXPECT_NEAR(desc[1], (r01 - r02) * (r01 - r02), 1e-9);
    EXPECT_NEAR(desc[2], r12, 1e-9);
    EXPECT_NEAR(desc[3], val01 * val02, 1e-9);

    // 4. Test evaluateAndDifferentiate()
    auto desc_and_derivs = trans.evaluateAndDifferentiate(triplet);
    // Check value part
    EXPECT_NEAR(desc_and_derivs.value[0], r01 + r02, 1e-9);
    EXPECT_NEAR(desc_and_derivs.value[1], (r01 - r02) * (r01 - r02), 1e-9);
    EXPECT_NEAR(desc_and_derivs.value[2], r12, 1e-9);
    EXPECT_NEAR(desc_and_derivs.value[3], val01 * val02, 1e-9);

    // Check derivatives part
    const auto& derivs = desc_and_derivs.derivatives;
    EXPECT_EQ(derivs.size(), 3); // 3 separations in a triplet

    // Derivatives wrt r01
    EXPECT_NEAR(derivs[0][0], 1.0, 1e-9);
    EXPECT_NEAR(derivs[0][1], 2.0 * (r01 - r02), 1e-9);
    EXPECT_NEAR(derivs[0][2], 0.0, 1e-9);
    EXPECT_NEAR(derivs[0][3], deriv01 * val02, 1e-9);

    // Derivatives wrt r02
    EXPECT_NEAR(derivs[1][0], 1.0, 1e-9);
    EXPECT_NEAR(derivs[1][1], 2.0 * (r02 - r01), 1e-9);
    EXPECT_NEAR(derivs[1][2], 0.0, 1e-9);
    EXPECT_NEAR(derivs[1][3], val01 * deriv02, 1e-9);

    // Derivatives wrt r12
    EXPECT_NEAR(derivs[2][0], 0.0, 1e-9);
    EXPECT_NEAR(derivs[2][1], 0.0, 1e-9);
    EXPECT_NEAR(derivs[2][2], 1.0, 1e-9);
    EXPECT_NEAR(derivs[2][3], 0.0, 1e-9);
}
