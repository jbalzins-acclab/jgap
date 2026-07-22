#include "HermiteCubicSpline.hpp"

#include <algorithm>

#include "io/log/CurrentLogger.hpp"

namespace jgap {
    HermiteCubicSpline::HermiteCubicSpline(const Grid<1>& table) { init(table); }

    InterpolationResults<1> HermiteCubicSpline::interpolate(std::array<Real, 1> pos) const {
        const Real r = pos[0];
        const Real origin = table.origin[0];
        const Real spacing = table.spacing[0];
        const size_t n = table.sizes[0];

        if (r < origin) return {table({0}), {0.0}};
        if (r > getCutoff()[0]) return {0.0, {0.0}};

        size_t i = findInterval(r);
        Real dx = r - (origin + static_cast<Real>(i) * spacing);

        Real value = table({i}) + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
        Real derivative = b[i] + 2.0 * c[i] * dx + 3.0 * d[i] * dx * dx;

        return {value, {derivative}};
    }

    std::array<Real, 1> HermiteCubicSpline::getCutoff() const {
        return {table.origin[0] + static_cast<Real>(table.sizes[0] - 1) * table.spacing[0]};
    }

    void HermiteCubicSpline::init(const Grid<1>& table_in) {
        table = table_in;
        const size_t n = table.sizes[0];
        const Real spacing = table.spacing[0];
        const Real inverse_spacing = 1.0 / spacing;

        if (n < 2) {
            JGAP_LOG_ERROR("Spline table must have at least 2 points.", true);
        }

        b.resize(n); // Slopes at each knot point
        c.resize(n - 1);
        d.resize(n - 1);

        // Step 1: Compute the local slopes (b) at each point using finite differences
        // Since we have a uniform grid from Table<1>, h_left == h_right == spacing
        for (size_t i = 1; i < n - 1; ++i) {
            Real delta_left = (table({i}) - table({i - 1})) * inverse_spacing;
            Real delta_right = (table({i + 1}) - table({i})) * inverse_spacing;

            // Weighted average based on interval sizes (reduces to standard centered difference on uniform grids)
            b[i] = (delta_left + delta_right) / 2.0;
        }

        // Boundary conditions for slopes: One-sided finite differences
        b[0] = (table({1}) - table({0})) * inverse_spacing; // Left boundary slope
        b[n - 1] = (table({n - 1}) - table({n - 2})) * inverse_spacing; // Right boundary slope

        // Step 2: Compute c and d coefficients for each interval using Hermite formulas
        for (size_t i = 0; i < n - 1; ++i) {
            Real delta = (table({i + 1}) - table({i})) * inverse_spacing;

            // Mathematical transformation solving for continuous values and slopes at the boundaries
            c[i] = (3.0 * delta - 2.0 * b[i] - b[i + 1]) * inverse_spacing;
            d[i] = (b[i] + b[i + 1] - 2.0 * delta) * inverse_spacing * inverse_spacing;
        }
    }

    size_t HermiteCubicSpline::findInterval(Real r) const {
        const Real origin = table.origin[0];
        const Real inverse_spacing = 1.0 / table.spacing[0];
        const size_t n = table.sizes[0];

        if (r < origin) return 0;
        if (r > getCutoff()[0]) return n - 2;

        size_t idx = static_cast<size_t>((r - origin) * inverse_spacing);
        return std::min(idx, n - 2);
    }
}
