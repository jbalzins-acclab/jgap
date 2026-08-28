#include "NaturalCubicSpline.hpp"

#include <algorithm>

#include "../io/log/CurrentLogger.hpp"

namespace jgap {
    NaturalCubicSpline::NaturalCubicSpline(const std::vector<Real>& r_vec_in, const std::vector<Real>& e_vec_in) {
        init(r_vec_in, e_vec_in);
    }

    InterpolationResults<1> NaturalCubicSpline::interpolate(std::array<Real, 1> pos) const {
        const Real r = pos[0];
        if (r < r_vec.front()) return {energies.front(), {0.0}};
        if (r > r_vec.back()) return {0.0, {0.0}};

        size_t i = findInterval(r);
        double dx = r - r_vec[i];

        Real value = energies[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
        Real derivative = b[i] + 2.0 * c[i] * dx + 3.0 * d[i] * dx * dx;

        return {value, {derivative}};
    }

    void NaturalCubicSpline::init(const std::vector<Real>& r, const std::vector<Real>& e) {
        if (r.size() != e.size() || r.size() < 2) {
            JGAP_LOG_ERROR("Spline reference vectors must be the same size and have at least 2 points.", true);
        }
        if (!std::ranges::is_sorted(r)) {
            JGAP_LOG_ERROR("Spline reference distances must be sorted", true);
        }

        r_vec = r;
        const size_t n = r.size();
        energies = e;
        b.resize(n - 1);
        c.resize(n);
        d.resize(n - 1);

        std::vector<double> h(n - 1);
        for (size_t i = 0; i < n - 1; ++i) h[i] = r[i + 1] - r[i];

        std::vector<double> alpha(n - 1);
        for (size_t i = 1; i < n - 1; i++) {
            alpha[i] =
                (3.0 / h[i]) * (energies[i + 1] - energies[i]) - (3.0 / h[i - 1]) * (energies[i] - energies[i - 1]);
        }

        std::vector<double> l(n), mu(n), z(n);
        l[0] = 1.0;
        mu[0] = z[0] = 0.0;

        for (size_t i = 1; i < n - 1; i++) {
            l[i] = 2.0 * (r_vec[i + 1] - r_vec[i - 1]) - h[i - 1] * mu[i - 1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
        }

        l[n - 1] = 1.0;
        z[n - 1] = c[n - 1] = 0.0;

        for (size_t j = n - 1; j-- > 0;) {
            c[j] = z[j] - mu[j] * c[j + 1];
            b[j] = (energies[j + 1] - energies[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0;
            d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
        }
    }

    size_t NaturalCubicSpline::findInterval(Real r) const {
        const auto it = std::ranges::upper_bound(r_vec, r);
        const size_t idx = std::max(size_t{}, static_cast<size_t>(it - r_vec.begin()) - 1);
        return std::min(idx, r_vec.size() - 2);
    }
}
