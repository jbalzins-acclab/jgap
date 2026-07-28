#include "CubicBSpline.hpp"

#include <algorithm>
#include <cmath>

namespace jgap {

    std::vector<Real> CubicBSpline::toSplineCoefficients(const std::vector<Real>& original, const Real spacing) {
        static constexpr std::array basis = {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0};

        const size_t n_coefficients = original.size() + 2;

        const Real inverse_spacing = 1.0 / spacing;
        const Real inverse_spacing_sq = inverse_spacing * inverse_spacing;
        std::vector<std::array<Real, 4>> bands(n_coefficients);
        bands[0] = {inverse_spacing_sq, -2.0_r * inverse_spacing_sq, inverse_spacing_sq, 0.0};
        bands[n_coefficients - 1] = {inverse_spacing_sq, -2.0_r * inverse_spacing_sq, inverse_spacing_sq, 0.0};

        for (size_t i = 1; i < n_coefficients - 1; i++) {
            bands[i][0] = basis[0];
            bands[i][1] = basis[1];
            bands[i][2] = basis[2];
            bands[i][3] = original[i - 1];
        }

        bands[0][1] /= bands[0][0];
        bands[0][2] /= bands[0][0];
        bands[0][3] /= bands[0][0];
        bands[0][0] = 1.0;
        bands[1][1] -= bands[1][0] * bands[0][1];
        bands[1][2] -= bands[1][0] * bands[0][2];
        bands[1][3] -= bands[1][0] * bands[0][3];
        bands[0][0] = 0;
        bands[1][2] /= bands[1][1];
        bands[1][3] /= bands[1][1];
        bands[1][1] = 1.0;

        for (size_t i = 2; i < n_coefficients - 1; i++) {
            bands[i][1] -= bands[i][0] * bands[i - 1][2];
            bands[i][3] -= bands[i][0] * bands[i - 1][3];
            bands[i][2] /= bands[i][1];
            bands[i][3] /= bands[i][1];
            bands[i][0] = 0.0;
            bands[i][1] = 1.0;
        }

        bands[n_coefficients - 1][1] -= bands[n_coefficients - 1][0] * bands[n_coefficients - 3][2];
        bands[n_coefficients - 1][3] -= bands[n_coefficients - 1][0] * bands[n_coefficients - 3][3];
        bands[n_coefficients - 1][2] -= bands[n_coefficients - 1][1] * bands[n_coefficients - 2][2];
        bands[n_coefficients - 1][3] -= bands[n_coefficients - 1][1] * bands[n_coefficients - 2][3];
        bands[n_coefficients - 1][3] /= bands[n_coefficients - 1][2];
        bands[n_coefficients - 1][2] = 1.0;

        std::vector<Real> coefficients(n_coefficients, 0.0);
        coefficients[n_coefficients - 1] = bands[n_coefficients - 1][3];
        for (size_t i = n_coefficients - 2; i > 0; i--) {
            coefficients[i] = bands[i][3] - bands[i][2] * coefficients[i + 1];
        }
        coefficients[0] = bands[0][3] - bands[0][1] * coefficients[1] - bands[0][2] * coefficients[2];

        return coefficients;
    }

    CubicBSpline::CubicBSpline(const Grid<1>& coefficients) : coefficients(coefficients) {}

    ValuePtr<Spline<1>> CubicBSpline::fit(const Grid<1>& values) {
        auto new_dims = values.sizes;
        new_dims[0] += 2;
        auto new_origin = values.origin;
        new_origin[0] -= values.spacing[0];

        Grid<1> coeff_table(new_dims, values.spacing, new_origin);
        coeff_table.data_flat = toSplineCoefficients(values.data_flat, values.spacing[0]);
        return CubicBSpline(coeff_table);
    }

    InterpolationResults<1> CubicBSpline::interpolate(std::array<Real, 1> pos) const {
        const Real r = pos[0];
        const Real h = coefficients.spacing[0];
        const Real data_origin = coefficients.origin[0] + h;
        const Real data_upper = data_origin + static_cast<Real>(coefficients.sizes[0] - 3) * h;
        if (r < (data_origin - 1e-9) || r > (data_upper + 1e-9)) {
            return {0.0, {0.0}};
        }

        const Real u = (r - data_origin) / h;
        // Maximum allowed starting index: coeff sizes - 4 ensures reading i..i+3 remains within valid memory [0,
        // sizes-1]
        const int imax = static_cast<int>(coefficients.sizes[0]) - 4;
        const int i = std::clamp(static_cast<int>(std::floor(u)), 0, imax);
        const Real t = u - i;

        const Real t2 = t * t;
        const Real t3 = t2 * t;

        Real w[4];
        w[0] = (1 - 3 * t + 3 * t2 - t3) / 6.0;
        w[1] = (4 - 6 * t2 + 3 * t3) / 6.0;
        w[2] = (1 + 3 * t + 3 * t2 - 3 * t3) / 6.0;
        w[3] = t3 / 6.0;

        Real dw[4];
        const Real dinv = 1.0 / h;
        dw[0] = (-3 + 6 * t - 3 * t2) * dinv / 6.0;
        dw[1] = (-12 * t + 9 * t2) * dinv / 6.0;
        dw[2] = (3 + 6 * t - 9 * t2) * dinv / 6.0;
        dw[3] = (3 * t2) * dinv / 6.0;

        Real value = 0.0;
        Real derivative = 0.0;
        for (int k = 0; k < 4; ++k) {
            value += coefficients.data_flat[i + k] * w[k];
            derivative += coefficients.data_flat[i + k] * dw[k];
        }

        return {value, {derivative}};
    }

    std::array<Real, 1> CubicBSpline::getCutoff() const {
        // The valid data range is from the original grid's origin (data_origin = coefficients.origin[0] + h)
        // to its last point at data_origin + (N - 1) * h, where N = coefficients.sizes[0] - 2.
        const Real data_points = static_cast<Real>(coefficients.sizes[0] - 2);
        const Real data_origin = coefficients.origin[0] + coefficients.spacing[0];
        return {data_origin + (data_points - 1) * coefficients.spacing[0]};
    }
}
