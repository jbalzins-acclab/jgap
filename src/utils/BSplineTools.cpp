#include "utils/BSplineTools.hpp"

namespace jgap {

    vector<double> BSplineTools::toSplineCoefficients(const vector<double> &original, const double spacing) {
        static constexpr array basis = {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0};

        const size_t nCoefficients = original.size() + 2;

        const double inverseSpacing = 1.0 / spacing;
        const double inverseSpacingSq = inverseSpacing * inverseSpacing;
        vector<array<double, 4>> bands(nCoefficients);
        bands[0] = {inverseSpacingSq, -2.0 * inverseSpacingSq, inverseSpacingSq, 0.0};
        bands[nCoefficients-1] = {inverseSpacingSq, -2.0 * inverseSpacingSq, inverseSpacingSq, 0.0};

        for (size_t i = 1; i < nCoefficients-1; i++) {
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

        for (size_t i = 2; i < nCoefficients-1; i++) {
            bands[i][1] -= bands[i][0] * bands[i-1][2];
            bands[i][3] -= bands[i][0] * bands[i-1][3];
            bands[i][2] /= bands[i][1];
            bands[i][3] /= bands[i][1];
            bands[i][0] = 0.0;
            bands[i][1] = 1.0;
        }

        bands[nCoefficients-1][1] -= bands[nCoefficients-1][0] * bands[nCoefficients-3][2];
        bands[nCoefficients-1][3] -= bands[nCoefficients-1][0] * bands[nCoefficients-3][3];
        bands[nCoefficients-1][2] -= bands[nCoefficients-1][1] * bands[nCoefficients-2][2];
        bands[nCoefficients-1][3] -= bands[nCoefficients-1][1] * bands[nCoefficients-2][3];
        bands[nCoefficients-1][3] /= bands[nCoefficients-1][2];
        bands[nCoefficients-1][2] = 1.0;

        vector coefficients(nCoefficients, 0.0);
        coefficients[nCoefficients-1] = bands[nCoefficients-1][3];
        for (size_t i = nCoefficients-2; i > 0; i--) {
            coefficients[i] = bands[i][3] - bands[i][2] * coefficients[i+1];
        }
        coefficients[0] = bands[0][3] - bands[0][1] * coefficients[1] - bands[0][2] * coefficients[2];

        return coefficients;
    }

    Grid1d BSplineTools::toSplineCoefficients(const Grid1d &original) {
        return Grid1d(
            toSplineCoefficients(original.data, original.spacing),
            original.spacing,
            original.origin - original.spacing
            );
    }

    Grid3d BSplineTools::toSplineCoefficients(const Grid3d &original) {
        Grid3d coeffGrid(
            original.originR() - original.spacingR(),
            original.nR + 2,
            original.spacingR(),
            original.originAngular() - original.spacingAngular(),
            original.nAngular + 2,
            original.spacingAngular()
            );

        for (size_t i = 0; i < original.nR; i++) {
            for (size_t j = 0; j < original.nR; j++) {
                const auto c = toSplineCoefficients(coeffGrid.slice_ij(i, j), original.spacingAngular());
                for (size_t k = 0; k < c.size(); k++) {
                    coeffGrid(i+1, j+1, k) = c[k];
                }
            }
        }

        for (size_t i = 0; i < original.nR; i++) {
            for (size_t k = 0; k < original.nAngular; k++) {
                const auto c = toSplineCoefficients(coeffGrid.slice_ik(i, k), original.spacingR());
                for (size_t j = 0; j < c.size(); j++) {
                    coeffGrid(i+1, j, k) = c[j];
                }
            }
        }
        for (size_t j = 0; j < original.nR; j++) {
            for (size_t k = 0; k < original.nAngular; k++) {
                const auto c = toSplineCoefficients(coeffGrid.slice_jk(j, k), original.spacingR());
                for (size_t i = 0; i < c.size(); i++) {
                    coeffGrid(i, j, k) = c[i];
                }
            }
        }

        return coeffGrid;
    }

    BSplineTools::InterpolationResults<double> BSplineTools::interpolate(const Grid1d &table, const double pos) {
        const auto &coeffs = table.data;
        const double h = table.spacing;
        const double origin = table.origin;

        // Find normalized coordinate
        const double x = (pos - origin) / h;
        const int i = static_cast<int>(floor(x));
        const double t = x - i;

        // Precompute basis weights
        double w[4];
        const double t2 = t * t;
        const double t3 = t2 * t;
        w[0] = (1 - 3*t + 3*t2 - t3) / 6.0;
        w[1] = (4 - 6*t2 + 3*t3) / 6.0;
        w[2] = (1 + 3*t + 3*t2 - 3*t3) / 6.0;
        w[3] = t3 / 6.0;

        // Interpolate value
        double value = 0.0;
        for (int k = 0; k < 4; ++k)
            value += coeffs[i + k - 1] * w[k];

        // Optional derivative (for gradient info)
        double dw[4];
        dw[0] = (-3 + 6*t - 3*t2) / 6.0;
        dw[1] = (-12*t + 9*t2) / 6.0;
        dw[2] = (3 + 6*t - 9*t2) / 6.0;
        dw[3] = (3*t2) / 6.0;

        double derivative = 0.0;
        for (int k = 0; k < 4; ++k)
            derivative += coeffs[i + k - 1] * dw[k];
        derivative /= h;

        return { value, derivative };
    }

    BSplineTools::InterpolationResults<Vector3> BSplineTools::interpolate(const Grid3d &table, const Vector3 &pos) {
        const double hx = table.spacingR();
        const double hy = table.spacingR();       // or spacingTheta if different
        const double hz = table.spacingAngular();

        // Normalized coordinates
        const double x = (pos.x - table.originR()) / hx;
        const double y = (pos.y - table.originAngular()) / hy;
        const double z = (pos.z - table.originAngular()) / hz;

        const int ix = static_cast<int>(floor(x));
        const int iy = static_cast<int>(floor(y));
        const int iz = static_cast<int>(floor(z));

        const double tx = x - ix;
        const double ty = y - iy;
        const double tz = z - iz;

        // Compute cubic B-spline basis and derivatives
        auto computeWeights = [](double t, double w[4], double dw[4]) {
            const double t2 = t * t;
            const double t3 = t2 * t;

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

        double wx[4], wy[4], wz[4];
        double dwx[4], dwy[4], dwz[4];
        computeWeights(tx, wx, dwx);
        computeWeights(ty, wy, dwy);
        computeWeights(tz, wz, dwz);

        // Accumulators
        double value = 0.0;
        Vector3 grad(0.0, 0.0, 0.0);

        // Tensor-product evaluation
        for (int a = 0; a < 4; ++a) {
            for (int b = 0; b < 4; ++b) {
                for (int c = 0; c < 4; ++c) {
                    const double coeff = table(ix + a - 1, iy + b - 1, iz + c - 1);

                    const double wabc = wx[a] * wy[b] * wz[c];
                    value += coeff * wabc;

                    grad.x += coeff * dwx[a] * wy[b] * wz[c];
                    grad.y += coeff * wx[a] * dwy[b] * wz[c];
                    grad.z += coeff * wx[a] * wy[b] * dwz[c];
                }
            }
        }

        // Scale derivatives by grid spacing
        grad.x /= hx;
        grad.y /= hy;
        grad.z /= hz;

        return { value, grad };
    }
}
