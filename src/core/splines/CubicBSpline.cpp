#include "CubicBSpline.hpp"

#include <cmath>

namespace jgap {

    std::vector<Real> CubicBSpline::toSplineCoefficients(const std::vector<Real> &original, const Real spacing) {
        static constexpr std::array basis = {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0};

        const size_t n_coefficients = original.size() + 2;

        const Real inverse_spacing = 1.0 / spacing;
        const Real inverse_spacing_sq = inverse_spacing * inverse_spacing;
        std::vector<std::array<Real, 4>> bands(n_coefficients);
        bands[0] = {inverse_spacing_sq, -2.0 * inverse_spacing_sq, inverse_spacing_sq, 0.0};
        bands[n_coefficients-1] = {inverse_spacing_sq, -2.0 * inverse_spacing_sq, inverse_spacing_sq, 0.0};

        for (size_t i = 1; i < n_coefficients-1; i++) {
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

        for (size_t i = 2; i < n_coefficients-1; i++) {
            bands[i][1] -= bands[i][0] * bands[i-1][2];
            bands[i][3] -= bands[i][0] * bands[i-1][3];
            bands[i][2] /= bands[i][1];
            bands[i][3] /= bands[i][1];
            bands[i][0] = 0.0;
            bands[i][1] = 1.0;
        }

        bands[n_coefficients-1][1] -= bands[n_coefficients-1][0] * bands[n_coefficients-3][2];
        bands[n_coefficients-1][3] -= bands[n_coefficients-1][0] * bands[n_coefficients-3][3];
        bands[n_coefficients-1][2] -= bands[n_coefficients-1][1] * bands[n_coefficients-2][2];
        bands[n_coefficients-1][3] -= bands[n_coefficients-1][1] * bands[n_coefficients-2][3];
        bands[n_coefficients-1][3] /= bands[n_coefficients-1][2];
        bands[n_coefficients-1][2] = 1.0;

        std::vector<Real> coefficients(n_coefficients, 0.0);
        coefficients[n_coefficients-1] = bands[n_coefficients-1][3];
        for (size_t i = n_coefficients-2; i > 0; i--) {
            coefficients[i] = bands[i][3] - bands[i][2] * coefficients[i+1];
        }
        coefficients[0] = bands[0][3] - bands[0][1] * coefficients[1] - bands[0][2] * coefficients[2];

        return coefficients;
    }

    CubicBSpline::CubicBSpline(const Table<1>& coefficients) : coefficients(coefficients) {}

    CubicBSpline CubicBSpline::fit(const Table<1>& values) {
        auto new_dims = values.dims;
        new_dims[0] += 2;
        auto new_origin = values.origin;
        new_origin[0] -= values.spacing[0];

        Table<1> coeff_table(new_dims, values.spacing, new_origin);
        coeff_table.data_flat = toSplineCoefficients(values.data_flat, values.spacing[0]);
        return CubicBSpline(coeff_table);
    }

    InterpolationResults<1> CubicBSpline::interpolate(std::array<Real, 1> pos) const {
        const Real r = pos[0];
        const auto &coeffs = coefficients.data_flat;
        const Real h = coefficients.spacing[0];
        const Real origin = coefficients.origin[0];

        // Find normalized coordinate
        const Real x = (r - origin) / h;
        const int i = static_cast<int>(std::floor(x));
        const Real t = x - i;

        // Precompute basis weights
        Real w[4];
        const Real t2 = t * t;
        const Real t3 = t2 * t;
        w[0] = (1 - 3*t + 3*t2 - t3) / 6.0;
        w[1] = (4 - 6*t2 + 3*t3) / 6.0;
        w[2] = (1 + 3*t + 3*t2 - 3*t3) / 6.0;
        w[3] = t3 / 6.0;

        // Interpolate value
        Real value = 0.0;
        for (int k = 0; k < 4; ++k)
            value += coeffs[i + k - 1] * w[k];

        // Optional derivative (for gradient info)
        Real dw[4];
        dw[0] = (-3 + 6*t - 3*t2) / 6.0;
        dw[1] = (-12*t + 9*t2) / 6.0;
        dw[2] = (3 + 6*t - 9*t2) / 6.0;
        dw[3] = (3*t2) / 6.0;

        Real derivative = 0.0;
        for (int k = 0; k < 4; ++k)
            derivative += coeffs[i + k - 1] * dw[k];
        derivative /= h;

        return { value, {derivative} };
    }

    std::array<Real, 1> CubicBSpline::getCutoff() const {
        return {coefficients.origin[0] + static_cast<Real>(coefficients.dims[0] - 2) * coefficients.spacing[0]};
    }
}