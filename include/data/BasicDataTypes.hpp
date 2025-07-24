#ifndef BASIC_DATA_TYPES_HPP
#define BASIC_DATA_TYPES_HPP

#include <algorithm>
#include <string>
#include <optional>
#include <vector>
#include <array>
#include <cmath>
#include <format>

#include "io/log/CurrentLogger.hpp"

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
        Vector3& operator+=(const Vector3& other) {
            x += other.x;
            y += other.y;
            z += other.z;
            return *this;
        }
        Vector3 operator-(const Vector3& other) const {
            return Vector3{x - other.x, y - other.y, z - other.z};
        };
        Vector3& operator-=(const Vector3& other) {
            x -= other.x;
            y -= other.y;
            z -= other.z;
            return *this;
        }
        Vector3 operator*(const double scalar) const {
            return Vector3{x * scalar, y * scalar, z * scalar};
        };
        Vector3& operator*=(const double scalar) {
            x *= scalar;
            y *= scalar;
            z *= scalar;
            return *this;
        }
        Vector3 operator/(const double scalar) const {
            return Vector3{x / scalar, y / scalar, z / scalar};
        }
        Vector3& operator/=(const double scalar) {
            x /= scalar;
            y /= scalar;
            z /= scalar;
            return *this;
        }
        double dot(const Vector3& other) const {
            return x * other.x + y * other.y + z * other.z;
        };
        Vector3 cross(const Vector3& other) const {
            return Vector3{
                y * other.z - z * other.y,
                z * other.x - x * other.z,
                x * other.y - y * other.x
            };
        }
        double square() const {
            return x * x + y * y + z * z;
        }
        double len() const {
            return sqrt(x * x + y * y + z * z);
        };
        double project(const Vector3& other) const {
            return dot(other) / other.len();
        }
        double aproject(const Vector3& other) const {
            return sqrt(len() * len() - project(other) * project(other));
        }
        Vector3 normalize() const {
            return *this * (1.0 / len());
        };
        double min() const {
            const double t = abs(x) < abs(y) ? x : y;
            return abs(t) < abs(z) ? abs(t) : abs(z);
        }
        string toString() const {
            return to_string(x) + ", " + to_string(y) + ", " + to_string(z);
        }
        double aproject(const Vector3& u, const Vector3& v) const {
            Vector3 _cross = u.cross(v);
            if (_cross.len() == 0.0) return this->aproject(u);
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

    struct PotentialPrediction {
        optional<double> energy;
        optional<vector<Vector3>> forces;
        optional<array<Vector3, 3>> virials;

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
            optional<array<Vector3, 3>> _virials;
            if (virials.has_value() && other.virials.has_value()) {
                _virials = virials;
                (*_virials)[0] += other.virials.value()[0];
                (*_virials)[1] += other.virials.value()[1];
                (*_virials)[2] += other.virials.value()[2];
            } else if (virials.has_value()) {
                _virials = virials.value();
            } else if (other.virials.has_value()) {
                _virials = other.virials.value();
            }
            return PotentialPrediction{
                _energy, _forces, _virials
            };
        }
    };

    using NeighboursData = vector<NeighbourData>;

    struct AtomicStructure {
        optional<string> configType;
        array<Vector3, 3> lattice;
        vector<Vector3> positions;
        vector<Species> species;

        optional<vector<NeighboursData>> neighbours;

        optional<double> energy;
        optional<vector<Vector3>> forces;
        optional<array<Vector3, 3>> virials;

        optional<double> energySigmaInverse;
        optional<vector<Vector3>> forceSigmasInverse;
        optional<array<Vector3, 3>> virialSigmasInverse;

        struct AtomProxy {
            size_t index;
            AtomicStructure* structure;

            Vector3& position() const { return structure->positions[index]; }
            Species& species() const { return structure->species[index]; }

            Vector3& force() const {
                if (!structure->forces)
                    throw std::runtime_error("Forces not set");
                return (*structure->forces)[index];
            }
            Vector3& forceSigmasInverse() const {
                if (!structure->forceSigmasInverse)
                    throw std::runtime_error("Forces not set");
                return (*structure->forceSigmasInverse)[index];
            }
            NeighboursData& neighbours() const {
                if (!structure->neighbours)
                    throw std::runtime_error("Neighbours not set");
                return (*structure->neighbours)[index];
            }
        };

        struct ConstAtomProxy {
            size_t index;
            const AtomicStructure* structure;

            Vector3 position() const { return structure->positions[index]; }
            Species species() const { return structure->species[index]; }

            Vector3 force() const {
                if (!structure->forces)
                    throw std::runtime_error("Forces not set");
                return (*structure->forces)[index];
            }
            Vector3 forceSigmasInverse() const {
                if (!structure->forceSigmasInverse)
                    throw std::runtime_error("Forces not set");
                return (*structure->forceSigmasInverse)[index];
            }
            NeighboursData neighbours() const {
                if (!structure->neighbours)
                    throw std::runtime_error("Neighbours not set");
                return (*structure->neighbours)[index];
            }
        };

        struct Iterator {
            AtomicStructure* structure;
            size_t index;

            Iterator& operator++() { ++index; return *this; }
            bool operator!=(const Iterator& other) const { return index != other.index; }
            AtomProxy operator*() const { return AtomProxy{index, structure}; }
        };

        struct ConstIterator {
            const AtomicStructure* structure;
            size_t index;

            ConstIterator& operator++() { ++index; return *this; }
            bool operator!=(const ConstIterator& other) const { return index != other.index; }
            ConstAtomProxy operator*() const { return ConstAtomProxy{index, structure}; }
        };

        Iterator begin() { return Iterator{this, 0}; }
        Iterator end() { return Iterator{this, positions.size()}; }

        ConstIterator begin() const { return ConstIterator{this, 0}; }
        ConstIterator end() const { return ConstIterator{this, positions.size()}; }

        AtomProxy operator[](const size_t i) { return AtomProxy{i, this}; }
        ConstAtomProxy operator[](const size_t i) const { return ConstAtomProxy{i, this}; }

        void setEnergyData(const PotentialPrediction& prediction);
        void adjust(const PotentialPrediction& prediction, bool subtract, bool setEmpty);
        AtomicStructure repeat(size_t a, size_t b, size_t c);

        double volume() const;
        size_t size() const { return positions.size(); }
    };

    struct Covariance {
        double total;
        vector<Vector3> forces;
        array<Vector3, 3> virials;
    };
}

#endif
