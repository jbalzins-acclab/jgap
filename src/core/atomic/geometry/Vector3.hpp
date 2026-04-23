#ifndef BASIC_DATA_TYPES_HPP
#define BASIC_DATA_TYPES_HPP

#include <string>
#include <optional>
#include <vector>
#include <array>
#include <cmath>
#include <map>

#include "../../Real.hpp"

namespace jgap {

    struct alignas(sizeof(Real) * 4) Vector3 {
    private:
        Real padding_ = 0;

    public:
        Vector3() : x(0.0), y(0.0), z(0.0) {}
        Vector3(Real x, Real y, Real z) : x(x), y(y), z(z) {}
        Vector3(const Vector3& other) = default;

        Real x, y, z;

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
        }
        Vector3& operator-=(const Vector3& other) {
            x -= other.x;
            y -= other.y;
            z -= other.z;
            return *this;
        }
        Vector3 operator*(const Real scalar) const {
            return Vector3{x * scalar, y * scalar, z * scalar};
        }
        Vector3 componentMul(const Vector3 other) const {
            return Vector3{x * other.x, y * other.y, z * other.z};
        };
        Vector3& operator*=(const Real scalar) {
            x *= scalar;
            y *= scalar;
            z *= scalar;
            return *this;
        }
        Vector3 operator/(const Real scalar) const {
            return Vector3{x / scalar, y / scalar, z / scalar};
        }
        Vector3& operator/=(const Real scalar) {
            x /= scalar;
            y /= scalar;
            z /= scalar;
            return *this;
        }
        Real dot(const Vector3& other) const {
            return x * other.x + y * other.y + z * other.z;
        };
        Vector3 cross(const Vector3& other) const {
            return Vector3{
                y * other.z - z * other.y,
                z * other.x - x * other.z,
                x * other.y - y * other.x
            };
        }
        Real square() const {
            return x * x + y * y + z * z;
        }
        Real norm() const {
            return sqrt(x * x + y * y + z * z);
        };
        Real project(const Vector3& other) const {
            return dot(other) / other.norm();
        }
        Real aproject(const Vector3& other) const {
            return sqrt(norm() * norm() - project(other) * project(other));
        }
        Vector3 normalize() const {
            return *this * (1.0 / norm());
        }
        Real min() const {
            const Real t = abs(x) < abs(y) ? x : y;
            return abs(t) < abs(z) ? abs(t) : abs(z);
        }
        Real aproject(const Vector3& u, const Vector3& v) const {
            Vector3 _cross = u.cross(v);
            if (_cross.norm() == 0.0) return this->aproject(u);
            return abs(this->project(_cross));
        }
        bool operator==(const Vector3& other) const {
            return x == other.x && y == other.y && z == other.z;
        }
        std::string toString() const {
            return std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z);
        }
    };
}

#endif
