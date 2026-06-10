#include "CubicBSpline3D.hpp"
#include "CubicBSpline.hpp"

namespace jgap {

    CubicBSpline3D::CubicBSpline3D(const Table<3>& coefficients) : coefficients(coefficients) {}

    CubicBSpline3D CubicBSpline3D::fit(const Table<3>& values) {
        auto new_dims = values.dims;
        new_dims[0] += 2;
        new_dims[1] += 2;
        new_dims[2] += 2;

        auto new_origin = values.origin;
        new_origin[0] -= values.spacing[0];
        new_origin[1] -= values.spacing[1];
        new_origin[2] -= values.spacing[2];

        Table<3> coeff_grid(new_dims, values.spacing, new_origin);

        for (size_t i = 0; i < values.dims[0]; i++) {
            for (size_t j = 0; j < values.dims[1]; j++) {
                std::vector<Real> slice(values.dims[2]);
                for(size_t k = 0; k < values.dims[2]; ++k) {
                    slice[k] = values({i, j, k});
                }
                const auto c = CubicBSpline::toSplineCoefficients(slice, values.spacing[2]);
                for (size_t k = 0; k < c.size(); k++) {
                    coeff_grid({i + 1, j + 1, k}) = c[k];
                }
            }
        }

        for (size_t i = 0; i < values.dims[0]; i++) {
            for (size_t k = 0; k < new_dims[2]; k++) {
                std::vector<Real> slice(values.dims[1]);
                for(size_t j = 0; j < values.dims[1]; ++j) {
                    slice[j] = coeff_grid({i + 1, j + 1, k});
                }
                const auto c = CubicBSpline::toSplineCoefficients(slice, values.spacing[1]);
                for (size_t j = 0; j < c.size(); j++) {
                    coeff_grid({i + 1, j, k}) = c[j];
                }
            }
        }
        for (size_t j = 0; j < new_dims[1]; j++) {
            for (size_t k = 0; k < new_dims[2]; k++) {
                std::vector<Real> slice(values.dims[0]);
                for(size_t i = 0; i < values.dims[0]; ++i) {
                    slice[i] = coeff_grid({i + 1, j, k});
                }
                const auto c = CubicBSpline::toSplineCoefficients(slice, values.spacing[0]);
                for (size_t i = 0; i < c.size(); i++) {
                    coeff_grid({i, j, k}) = c[i];
                }
            }
        }

        return CubicBSpline3D(coeff_grid);
    }

    InterpolationResults<3> CubicBSpline3D::interpolate(std::array<Real, 3> pos) const {
        const Real hx = coefficients.spacing[0];
        const Real hy = coefficients.spacing[1];
        const Real hz = coefficients.spacing[2];

        // Normalized coordinates
        const Real x = (pos[0] - coefficients.origin[0]) / hx;
        const Real y = (pos[1] - coefficients.origin[1]) / hy;
        const Real z = (pos[2] - coefficients.origin[2]) / hz;

        const size_t ix = static_cast<int>(floor(x));
        const size_t iy = static_cast<int>(floor(y));
        const size_t iz = static_cast<int>(floor(z));

        const Real tx = x - ix;
        const Real ty = y - iy;
        const Real tz = z - iz;

        // Compute cubic B-spline basis and derivatives
        auto compute_weights = [](Real t, Real w[4], Real dw[4]) {
            const Real t2 = t * t;
            const Real t3 = t2 * t;

            // Basis (B3)
            w[0] = (1.0 - 3.0*t + 3.0*t2 - t3) / 6.0;
            w[1] = (4.0 - 6.0*t2 + 3.0*t3) / 6.0;
            w[2] = (1.0 + 3.0*t + 3.0*t2 - 3.0*t3) / 6.0;
            w[3] = t3 / 6.0;

            // Derivative of basis (dB3/dt)
            dw[0] = (-3.0 + 6.0*t - 3.0*t2) / 6.0;
            dw[1] = (-12.0*t + 9.0*t2) / 6.0;
            dw[2] = (3.0 + 6.0*t - 9.0*t2) / 6.0;
            dw[3] = (3.0*t2) / 6.0;
        };

        Real wx[4], wy[4], wz[4];
        Real dwx[4], dwy[4], dwz[4];
        compute_weights(tx, wx, dwx);
        compute_weights(ty, wy, dwy);
        compute_weights(tz, wz, dwz);

        // Accumulators
        Real value = 0.0;
        std::array<Real, 3> grad = {0.0, 0.0, 0.0};

        // Tensor-product evaluation
        for (size_t a = 0; a < 4; ++a) {
            for (size_t b = 0; b < 4; ++b) {
                for (size_t c = 0; c < 4; ++c) {
                    const Real coeff = coefficients({(ix + a - 1), (iy + b - 1), (iz + c - 1)});

                    const Real wabc = wx[a] * wy[b] * wz[c];
                    value += coeff * wabc;

                    grad[0] += coeff * dwx[a] * wy[b] * wz[c];
                    grad[1] += coeff * wx[a] * dwy[b] * wz[c];
                    grad[2] += coeff * wx[a] * wy[b] * dwz[c];
                }
            }
        }

        // Scale derivatives by grid spacing
        grad[0] /= hx;
        grad[1] /= hy;
        grad[2] /= hz;

        return { value, grad };
    }

    std::array<Real, 3> CubicBSpline3D::getCutoff() const {
        return {coefficients.origin[0] + static_cast<Real>(coefficients.dims[0] - 2) * coefficients.spacing[0],
                coefficients.origin[1] + static_cast<Real>(coefficients.dims[1] - 2) * coefficients.spacing[1],
                coefficients.origin[2] + static_cast<Real>(coefficients.dims[2] - 2) * coefficients.spacing[2]};
    }
}