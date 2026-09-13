#include <gtest/gtest.h>
#include "jgap/core/atomic/species/Species.hpp"
#include "jgap/core/atomic/species/composition/Species2Sorted.hpp"
#include "jgap/core/potentials/coulomb/ScreenedCoulombPotential.hpp"
#include "jgap/core/splines/CubicBSpline.hpp"
#include "jgap/core/splines/CubicBSpline3D.hpp"
#include "jgap/core/splines/Grid.hpp"
#include "jgap/core/splines/HermiteCubicSpline.hpp"
#include "jgap/core/splines/NaturalCubicSpline.hpp"

#include <cmath>

using namespace jgap;

TEST(SplineTest, OneDimensionalSplines) {
    Species si("Si");
    ScreenedCoulombPotential sc({Species2Sorted(si, si)});

    std::vector<double> r_vec;
    std::vector<double> e_vec;
    std::vector<double> de_vec;

    double cutoff = 4.0;
    double spacing = 0.001;
    for (double r = 0.1; r <= cutoff + 1e-9; r += spacing) {
        r_vec.push_back(r);
        auto ed = sc.energyAndDerivative(Species2Sorted(si, si), r);
        e_vec.push_back(ed[0]);
        de_vec.push_back(ed[1]);
    }

    Grid<1> grid({e_vec.size()}, {spacing}, {0.1}, e_vec);

    NaturalCubicSpline natural_cubic_spline(r_vec, e_vec);
    HermiteCubicSpline hermite_cubic_spline(grid);
    auto b_spline_ptr = CubicBSpline::fit(grid);
    CubicBSpline& b_spline = *dynamic_cast<CubicBSpline*>(b_spline_ptr.get());

    for (size_t i = 0; i < r_vec.size(); ++i) {
        double r = r_vec[i];
        double sc_e = e_vec[i];
        double sc_de = de_vec[i];

        auto natural_res = natural_cubic_spline.interpolate({r});
        auto hermite_res = hermite_cubic_spline.interpolate({r});
        auto b_spline_res = b_spline.interpolate({r});

        EXPECT_NEAR(sc_e, natural_res.value, 1e-9);
        EXPECT_NEAR(sc_e, hermite_res.value, 1e-9);
        EXPECT_NEAR(sc_e, b_spline_res.value, 1e-9);

        double de_tolerance = 1e-9 + std::max(2e-5, std::abs(sc_de) * 1.2e-2);
        EXPECT_NEAR(sc_de, natural_res.gradient[0], de_tolerance);
        EXPECT_NEAR(sc_de, hermite_res.gradient[0], de_tolerance);
        EXPECT_NEAR(sc_de, b_spline_res.gradient[0], de_tolerance);
    }

    for (double r = 0.1005; r < cutoff; r += spacing) {
        auto ed = sc.energyAndDerivative(Species2Sorted(si, si), r);
        double sc_e = ed[0];
        double sc_de = ed[1];

        auto natural_res = natural_cubic_spline.interpolate({r});
        auto hermite_res = hermite_cubic_spline.interpolate({r});
        auto b_spline_res = b_spline.interpolate({r});

        double e_tolerance = 1e-9 + std::abs(sc_e) * 1e-3;
        EXPECT_NEAR(natural_res.value, sc_e, e_tolerance);
        EXPECT_NEAR(hermite_res.value, sc_e, e_tolerance);
        EXPECT_NEAR(b_spline_res.value, sc_e, e_tolerance);

        double de_tolerance = 1e-9 + std::max(2e-5, std::abs(sc_de) * 1e-2);
        EXPECT_NEAR(natural_res.gradient[0], sc_de, de_tolerance);
        EXPECT_NEAR(hermite_res.gradient[0], sc_de, de_tolerance);
        EXPECT_NEAR(b_spline_res.gradient[0], sc_de, de_tolerance);
    }
}

TEST(SplineTest, CubicBSplineCutoff) {
    std::vector<double> e_vec = {1.0, 2.0, 3.0, 2.0, 1.0, 0.0};
    Grid<1> grid({e_vec.size()}, {1.0}, {0.0}, e_vec);
    auto b_spline_ptr = CubicBSpline::fit(grid);
    CubicBSpline& b_spline = *dynamic_cast<CubicBSpline*>(b_spline_ptr.get());

    double spline_cutoff = b_spline.getCutoff()[0];
    EXPECT_NEAR(spline_cutoff, 5.0, 1e-9);

    // Test just inside the cutoff
    auto result_inside_cutoff = b_spline.interpolate({spline_cutoff - 1e-9});
    EXPECT_NEAR(result_inside_cutoff.value, 0.0, 2e-9);

    // Test at the cutoff
    // Since we relaxed the boundary check to +1e-9 to allow continuous extrapolation, 
    // the spline will evaluate its actual gradient here rather than an abrupt 0.0.
    auto result_at_cutoff = b_spline.interpolate({spline_cutoff});
    EXPECT_NEAR(result_at_cutoff.value, 0.0, 1e-9);
    EXPECT_NEAR(result_at_cutoff.gradient[0], result_inside_cutoff.gradient[0], 1e-5);

    // Test beyond the cutoff (outside the 1e-9 tolerance)
    auto result_beyond_cutoff = b_spline.interpolate({spline_cutoff + 1.0});
    EXPECT_NEAR(result_beyond_cutoff.value, 0.0, 1e-9);
    EXPECT_NEAR(result_beyond_cutoff.gradient[0], 0.0, 1e-9);
}

TEST(SplineTest, PreCutoffZero) {
    // Create data that is zero for the last few points
    std::vector<double> r_vec;
    std::vector<double> e_vec;
    double spacing = 0.5;
    for (double r = 0.0; r <= 5.0 + 1e-9; r += spacing) {
        r_vec.push_back(r);
        if (r < 3.0) {
            e_vec.push_back(std::cos(r * M_PI / 6.0)); // Function is zero at r=3
        } else {
            e_vec.push_back(0.0);
        }
    }

    // 1. NaturalCubicSpline
    NaturalCubicSpline natural_spline(r_vec, e_vec);
    EXPECT_NEAR(natural_spline.interpolate({3.0}).value, 0.0, 1e-9);
    EXPECT_NEAR(natural_spline.interpolate({3.5}).value, 0.0, 1e-9);
    EXPECT_NEAR(natural_spline.interpolate({4.0}).value, 0.0, 1e-9);

    // 2. Grid-based splines
    Grid<1> grid({e_vec.size()}, {spacing}, {0.0}, e_vec);

    // 2a. HermiteCubicSpline
    HermiteCubicSpline hermite_spline(grid);
    EXPECT_NEAR(hermite_spline.interpolate({3.0}).value, 0.0, 1e-9);
    EXPECT_NEAR(hermite_spline.interpolate({3.5}).value, 0.0, 1e-9);
    EXPECT_NEAR(hermite_spline.interpolate({4.0}).value, 0.0, 1e-9);

    // 2b. CubicBSpline
    auto b_spline_ptr = CubicBSpline::fit(grid);
    CubicBSpline& b_spline = *dynamic_cast<CubicBSpline*>(b_spline_ptr.get());
    EXPECT_NEAR(b_spline.interpolate({3.0}).value, 0.0, 1e-9);
    EXPECT_NEAR(b_spline.interpolate({3.5}).value, 0.0, 1e-9);
    EXPECT_NEAR(b_spline.interpolate({4.0}).value, 0.0, 1e-9);
}

TEST(SplineTest, ThreeDimensionalSpline) {
    auto mock_function = [](double x, double y, double z) { return std::sin(x) * std::cos(y) * std::exp(-z); };

    size_t n_points = 80;
    double spacing = 0.1;
    std::vector<double> values;
    for (size_t i = 0; i < n_points; ++i) {
        for (size_t j = 0; j < n_points; ++j) {
            for (size_t k = 0; k < n_points; ++k) {
                double x = 0.1 + i * spacing;
                double y = 0.1 + j * spacing;
                double z = 0.1 + k * spacing;
                values.push_back(mock_function(x, y, z));
            }
        }
    }

    Grid<3> grid({n_points, n_points, n_points}, {spacing, spacing, spacing}, {0.1, 0.1, 0.1}, values);
    CubicBSpline3D spline = CubicBSpline3D::fit(grid);

    double x = 0.15;
    double y = 0.25;
    double z = 0.35;
    double expected = mock_function(x, y, z);
    double interpolated = spline.interpolate({x, y, z}).value;

    EXPECT_NEAR(expected, interpolated, 1e-3);
}
