#ifndef BASIC_DATA_TYPES_HPP
#define BASIC_DATA_TYPES_HPP

#include <algorithm>
#include <string>
#include <optional>
#include <vector>
#include <array>
#include <cmath>
#include <format>

#include "io/log/Logger.hpp"

using namespace std;

namespace jgap {
    using Species = string;
    using OrderedSpeciesPair = pair<Species, Species>;

    class SpeciesPair {
    public:
        SpeciesPair(const Species& first, const Species& second)
            : pair_(minmax(first, second)) {}

        [[nodiscard]] const Species& first() const { return pair_.first; }
        [[nodiscard]] const Species& second() const { return pair_.second; }

        bool operator==(const SpeciesPair& other) const {
            return pair_ == other.pair_;
        }

        bool operator<(const SpeciesPair& other) const {
            return pair_ < other.pair_;
        }

        [[nodiscard]] string toString() const {
            return pair_.first + "," + pair_.second;
        }

        [[nodiscard]] bool contains(const Species& species) const {
            return pair_.first == species || pair_.second == species;
        }

    private:
        pair<Species, Species> pair_;
    };

    struct SpeciesTriplet {

        Species root;
        SpeciesPair nodes;

        bool operator==(const SpeciesTriplet& other) const {
            return root == other.root && nodes == other.nodes;
        }

        bool operator<(const SpeciesTriplet& other) const {
            if (root == other.root) {
                return nodes < other.nodes;
            }
            return root < other.root;
        }

        [[nodiscard]] string toString() const {
            return root + "," + nodes.toString();
        }
    };

    struct Vector3 {
        Vector3() = default;
        Vector3(double x, double y, double z) : x(x), y(y), z(z) {};
        Vector3(const Vector3& other) = default;

        double x, y, z;

        Vector3 operator+(const Vector3& other) const {
            return Vector3{x + other.x, y + other.y, z + other.z};
        }

        Vector3 operator-(const Vector3& other) const {
            return Vector3{x - other.x, y - other.y, z - other.z};
        };
        Vector3 operator*(const double scalar) const {
            return Vector3{x * scalar, y * scalar, z * scalar};
        };
        [[nodiscard]] double dot(const Vector3& other) const {
            return x * other.x + y * other.y + z * other.z;
        };
        [[nodiscard]] Vector3 cross(const Vector3& other) const {
            return Vector3{
                y * other.z - z * other.y,
                z * other.x - x * other.z,
                x * other.y - y * other.x
            };
        }
        [[nodiscard]] double square() const {
            return x * x + y * y + z * z;
        }
        [[nodiscard]] double norm() const {
            return sqrt(x * x + y * y + z * z);
        };
        [[nodiscard]] double project(const Vector3& other) const {
            return dot(other) / other.norm();
        }
        [[nodiscard]] double aproject(const Vector3& other) const {
            return sqrt(norm() * norm() - project(other) * project(other));
        }
        [[nodiscard]] Vector3 normalize() const {
            return *this * (1.0 / norm());
        };
        [[nodiscard]] double min() const {
            const double t = abs(x) < abs(y) ? x : y;
            return abs(t) < abs(z) ? abs(t) : abs(z);
        }
        [[nodiscard]] string toString() const {
            return to_string(x) + ", " + to_string(y) + ", " + to_string(z);
        }
        [[nodiscard]] double aproject(const Vector3& u, const Vector3& v) const {
            Vector3 _cross = u.cross(v);
            if (_cross.norm() == 0.0) return this->aproject(u);
            return abs(this->project(_cross));
        }
        bool operator==(const Vector3& other) const {
            return x == other.x && y == other.y && z == other.z;
        }
    };

    struct NeighbourData {
        size_t index;
        Vector3 offset;

        double distance;
    };

    struct AtomData {
        Vector3 position;
        Species species;

        // TODO: use KE+std(T) from relaxation in kernel?
        optional<Vector3> velocity; // why not
        optional<Vector3> force;
        optional<Vector3> forceSigmas;
        optional<vector<NeighbourData>> neighbours;
    };

    struct PotentialPrediction {
        optional<double> energy;
        optional<vector<Vector3>> forces;
        optional<array<array<double, 3>, 3>> virials;

        PotentialPrediction operator+(const PotentialPrediction& other) const {
            optional<double> _energy;
            if (energy.has_value() || other.energy.has_value()) {
                _energy = energy.value_or(0.0) + other.energy.value_or(0.0);
            }
            optional<vector<Vector3>> _forces;
            if (forces.has_value() && other.forces.has_value()) {
                _forces = vector<Vector3>(forces.value().size());
                for (size_t i = 0; i < forces.value().size(); i++) {
                    _forces->at(i) = forces.value()[i] + other.forces.value()[i];
                }
            } else if (forces.has_value()) {
                _forces = forces.value();
            } else if (other.forces.has_value()) {
                _forces = other.forces.value();
            }
            return PotentialPrediction{
                _energy, _forces,
            };
        }
    };

    struct AtomicStructure {
        array<Vector3, 3> latticeVectors;
        vector<AtomData> atoms;

        optional<string> configType;

        optional<double> energy;
        optional<double> energySigma;

        optional<array<array<double, 3>, 3>> virials;

        void set(const PotentialPrediction& prediction);
        void adjust(const PotentialPrediction& prediction, bool subtract, bool setEmpty);
        AtomicStructure repeat(size_t a, size_t b, size_t c);
    };

    struct Covariance {
        double total;
        vector<Vector3> derivatives;
    };
}

#endif
