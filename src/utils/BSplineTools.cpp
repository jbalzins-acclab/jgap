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

    BSplineTools::InterpolationResults<double> BSplineTools::interpolate(const Grid1d &table, double pos) {

        return Grid1d(coefficients, original.spacing, original.origin - original.spacing);
    }

    BSplineTools::InterpolationResults<Vector3> BSplineTools::interpolate(const Grid3d &table, Vector3 pos) {

    }
}
