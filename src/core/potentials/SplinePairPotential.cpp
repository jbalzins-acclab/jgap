#include "../../../include/core/potentials/SplinePairPotential.hpp"

#include <vector>
#include <stdexcept>
#include <iostream>
#include <algorithm>

#include "utils/Utils.hpp"

namespace jgap {
    SplinePairPotential::NaturalCubicSpline::NaturalCubicSpline(nlohmann::json params) {
        vector<double> distances{}, energies{};
        for (const double r: params["r"]) {
            distances.push_back(r);
        }
        for (const double E: params["E"]) {
            energies.push_back(E);
        }
        init(distances, energies);
    }

    SplinePairPotential::NaturalCubicSpline::NaturalCubicSpline(vector<double> r, vector<double> E) {
        init(r, E);
    }

    double SplinePairPotential::NaturalCubicSpline::evaluate(double r) const {
        if (r < _r.front() || r > _r.back())
            return 0.0;  // Optionally clamp or extrapolate

        size_t i = findInterval(r);
        double dx = r - _r[i];
        return _a[i] + _b[i] * dx + _c[i] * dx * dx + _d[i] * dx * dx * dx;
    }

    double SplinePairPotential::NaturalCubicSpline::derivative(double r) const {
        if (r < _r.front() || r > _r.back())
            return 0.0;

        size_t i = findInterval(r);
        double dx = r - _r[i];
        return _b[i] + 2.0 * _c[i] * dx + 3.0 * _d[i] * dx * dx;
    }

    nlohmann::json SplinePairPotential::NaturalCubicSpline::serialize() const {
        return nlohmann::json{
            {"r", _r},
            {"E", _a}
        };
    }

    void SplinePairPotential::NaturalCubicSpline::init(const vector<double> &r, const vector<double> &E) {
        if (r.size() != E.size() || r.size() < 2) {
            throw invalid_argument("Input vectors must be the same size and have at least 2 points.");
        }

        _r = r;
        const size_t n = r.size();
        _a = E;
        _b.resize(n - 1);
        _c.resize(n);
        _d.resize(n - 1);

        vector<double> h(n - 1);
        for (size_t i = 0; i < n - 1; ++i)
            h[i] = r[i + 1] - r[i];

        vector<double> alpha(n - 1);
        for (size_t i = 1; i < n - 1; ++i)
            alpha[i] = (3.0 / h[i]) * (_a[i + 1] - _a[i]) - (3.0 / h[i - 1]) * (_a[i] - _a[i - 1]);

        vector<double> l(n), mu(n), z(n);
        l[0] = 1.0;
        mu[0] = z[0] = 0.0;

        for (size_t i = 1; i < n - 1; ++i) {
            l[i] = 2.0 * (_r[i + 1] - _r[i - 1]) - h[i - 1] * mu[i - 1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
        }

        l[n - 1] = 1.0;
        z[n - 1] = _c[n - 1] = 0.0;

        for (int j = n - 2; j >= 0; --j) {
            _c[j] = z[j] - mu[j] * _c[j + 1];
            _b[j] = (_a[j + 1] - _a[j]) / h[j] - h[j] * (_c[j + 1] + 2.0 * _c[j]) / 3.0;
            _d[j] = (_c[j + 1] - _c[j]) / (3.0 * h[j]);
        }
    }

    size_t SplinePairPotential::NaturalCubicSpline::findInterval(double r) const {
        const auto it = ranges::upper_bound(_r, r);
        const size_t idx = max(static_cast<size_t>(0), static_cast<size_t>(it - _r.begin()) - 1);
        return min(idx, _r.size() - 2);
    }

    SplinePairPotential::SplinePairPotential(nlohmann::json params) {
        for (const auto &[speciesPairStr, pairParams]: params["pair_data"].items()) {
            Species species1 = split(speciesPairStr, ',')[0];
            Species species2 = split(speciesPairStr, ',')[1];
            _perSpeciesInterpolators[{species1, species2}] = make_shared<NaturalCubicSpline>(pairParams);
        }
    }

    SplinePairPotential::SplinePairPotential(map<SpeciesPair, pair<vector<double>, vector<double>>> points) {
        for (const auto &[speciesPair, pairParams]: points) {
            _perSpeciesInterpolators[speciesPair] = make_shared<NaturalCubicSpline>(
                pairParams.first, pairParams.second
                );
        }
    }

    PotentialPrediction SplinePairPotential::predict(const AtomicStructure &structure) {

        double energy = 0;
        vector forces(structure.size(), Vector3{0, 0, 0});
        array<Vector3, 3> virials{};

        for (size_t i = 0; i < structure.size(); i++) {

            auto atom1 = structure[i];

            for (const NeighbourData& neighbour: atom1.neighbours()) {

                auto atom2 = structure[neighbour.index];
                auto interpolator = _perSpeciesInterpolators[{atom1.species(), atom2.species()}];

                if (neighbour.index < i || neighbour.distance > interpolator->getCutoff()) continue;

                double dE = interpolator -> evaluate(neighbour.distance);
                double dE_dr = interpolator ->derivative(neighbour.distance);
                Vector3 r21 = atom1.position() - (atom2.position() + neighbour.offset);
                Vector3 f21 = r21.normalize() * -dE_dr;

                if (neighbour.index == i) {
                    dE /= 2.0;
                    f21 /= 2.0;
                } else {
                    forces[i] = forces[i] + f21;
                    forces[neighbour.index] = forces[neighbour.index] - f21;
                }

                energy += dE;
                virials[0] += f21 * r21.x;
                virials[1] += f21 * r21.y;
                virials[2] += f21 * r21.z;
            }
        }

        return PotentialPrediction{
            energy,
            forces,
            virials
        };
    }

    nlohmann::json SplinePairPotential::serialize() {
        nlohmann::json result{};

        for (const auto &[speciesPair, interpolator]: _perSpeciesInterpolators) {
            result[speciesPair.toString()] = interpolator->serialize();
        }

        return nlohmann::json{
            {"pair_data", result}
        };
    }

    double SplinePairPotential::getCutoff() {
        double cutoff = 0;
        for (const auto interpolator: _perSpeciesInterpolators | views::values) {
            cutoff = max(cutoff, interpolator->getCutoff());
        }
        return cutoff;
    }

    TabulationData SplinePairPotential::tabulate(const TabulationParams &params) {

        map<SpeciesPair, vector<double>> pairEnergies{};

        for (const auto &[speciesPair, interpolator]: _perSpeciesInterpolators) {
            vector<double> energies{};
            for (const double r: params.grid2b) {
                energies.push_back(interpolator->evaluate(r));
            }
            pairEnergies[speciesPair] = energies;
        }

        return TabulationData{.pairEnergies = pairEnergies};
    }
}
