#include "core/potentials/SplinePairPotential.hpp"

#include <vector>
#include <algorithm>

#include "io/log/CurrentLogger.hpp"
#include "utils/Utils.hpp"

namespace jgap {
    SplinePairPotential::NaturalCubicSpline::NaturalCubicSpline(DataNode params) {
        std::vector<double> distances{}, energies{};
        for (const auto& r_node : params["r"].asArray()) {
            distances.push_back(r_node.asDouble());
        }
        for (const auto& E_node : params["E"].asArray()) {
            energies.push_back(E_node.asDouble());
        }
        init(distances, energies);
    }

    SplinePairPotential::NaturalCubicSpline::NaturalCubicSpline(const std::vector<double>& r, const std::vector<double>& E) {
        init(r, E);
    }

    double SplinePairPotential::NaturalCubicSpline::evaluate(double r) const {
        if (r < r_.front()) return energies_.front();
        if (r > r_.back())  return 0.0; // cutoff behavior

        size_t i = findInterval(r);
        double dx = r - r_[i];
        return energies_[i] + b_[i] * dx + c_[i] * dx * dx + d_[i] * dx * dx * dx;
    }

    double SplinePairPotential::NaturalCubicSpline::derivative(double r) const {
        if (r < r_.front() || r > r_.back())
            return 0.0;

        size_t i = findInterval(r);
        double dx = r - r_[i];
        return b_[i] + 2.0 * c_[i] * dx + 3.0 * d_[i] * dx * dx;
    }

    DataNode SplinePairPotential::NaturalCubicSpline::serialize() const {
        return DataNode{
            {"r", r_},
            {"E", energies_}
        };
    }

    void SplinePairPotential::NaturalCubicSpline::init(const std::vector<double> &r, const std::vector<double> &E) {
        if (r.size() != E.size() || r.size() < 2) {
            JGAP_LOG_ERROR(
                "Spline reference vectors must be the same size and have at least 2 points.", true
                );
        }
        if (!ranges::is_sorted(r)) {
            JGAP_LOG_ERROR("Spline reference distances must be sorted", true);
        }

        r_ = r;
        const size_t n = r.size();
        energies_ = E;
        b_.resize(n - 1);
        c_.resize(n);
        d_.resize(n - 1);

        std::vector<double> h(n - 1);
        for (size_t i = 0; i < n - 1; ++i)
            h[i] = r[i + 1] - r[i];

        std::vector<double> alpha(n - 1);
        for (size_t i = 1; i < n - 1; ++i)
            alpha[i] = (3.0 / h[i]) * (energies_[i + 1] - energies_[i]) - (3.0 / h[i - 1]) * (energies_[i] - energies_[i - 1]);

        std::vector<double> l(n), mu(n), z(n);
        l[0] = 1.0;
        mu[0] = z[0] = 0.0;

        for (size_t i = 1; i < n - 1; ++i) {
            l[i] = 2.0 * (r_[i + 1] - r_[i - 1]) - h[i - 1] * mu[i - 1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
        }

        l[n - 1] = 1.0;
        z[n - 1] = c_[n - 1] = 0.0;

        for (int j = n - 2; j >= 0; --j) {
            c_[j] = z[j] - mu[j] * c_[j + 1];
            b_[j] = (energies_[j + 1] - energies_[j]) / h[j] - h[j] * (c_[j + 1] + 2.0 * c_[j]) / 3.0;
            d_[j] = (c_[j + 1] - c_[j]) / (3.0 * h[j]);
        }
    }

    size_t SplinePairPotential::NaturalCubicSpline::findInterval(double r) const {
        const auto it = ranges::upper_bound(r_, r);
        const size_t idx = std::max(static_cast<size_t>(0), static_cast<size_t>(it - r_.begin()) - 1);
        return std::min(idx, r_.size() - 2);
    }

    SplinePairPotential::SplinePairPotential(const DataNode& params) {
        for (const auto &[speciesPairStr, pairParams]: params["pair_data"].asObject()) {
            Species species1 = split(speciesPairStr, ',')[0];
            Species species2 = split(speciesPairStr, ',')[1];
            per_species_interpolators_[{species1, species2}] = std::make_shared<NaturalCubicSpline>(pairParams);
        }
    }

    SplinePairPotential::SplinePairPotential(const std::map<SpeciesPair, std::pair<std::vector<double>, std::vector<double>>>& points) {
        for (const auto &[speciesPair, pairParams]: points) {
            per_species_interpolators_[speciesPair] = std::make_shared<NaturalCubicSpline>(
                pairParams.first, pairParams.second
                );
        }
    }

    Predictions SplinePairPotential::predict(const AtomicStructure &structure) {

        double energy = 0;
        std::vector forces(structure.size(), Vector3{0, 0, 0});
        array<Vector3, 3> virials{};

        for (size_t i = 0; i < structure.size(); i++) {

            auto atom1 = structure[i];

            for (const NeighbourData& neighbour: atom1.neighbours()) {

                auto atom2 = structure[neighbour.index];
                auto interpolator = per_species_interpolators_[{atom1.species(), atom2.species()}];

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

        return Predictions{ energy, forces, virials };
    }

    DataNode SplinePairPotential::serialize() {
        DataNode result{};

        for (const auto &[speciesPair, interpolator]: per_species_interpolators_) {
            result[speciesPair.toString()] = interpolator->serialize();
        }

        return DataNode{
            {"pair_data", result}
        };
    }

    CutoffRanges SplinePairPotential::getCutoff() {
        double cutoff = 0;
        for (const auto& interpolator: per_species_interpolators_ | std::views::values) {
            cutoff = std::max(cutoff, interpolator->getCutoff());
        }
        return {.twoBody = cutoff};
    }

    void SplinePairPotential::tabulate(TabulationData &table) {
        for (const auto &[speciesPair, interpolator]: per_species_interpolators_) {
            for (const auto& it: table.getOrMake2bGrid(speciesPair)) {
                it.value += interpolator->evaluate(it.pos);
            }
        }
    }
}
