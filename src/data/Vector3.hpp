#ifndef BASIC_DATA_TYPES_HPP
#define BASIC_DATA_TYPES_HPP

#include <string>
#include <optional>
#include <vector>
#include <array>
#include <cmath>
#include <map>


namespace jgap {

    struct Vector3 {
        Vector3() : x(0.0), y(0.0), z(0.0) {}
        Vector3(double x, double y, double z) : x(x), y(y), z(z) {}
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
        }
        Vector3& operator-=(const Vector3& other) {
            x -= other.x;
            y -= other.y;
            z -= other.z;
            return *this;
        }
        Vector3 operator*(const double scalar) const {
            return Vector3{x * scalar, y * scalar, z * scalar};
        }
        Vector3 componentMul(const Vector3 other) const {
            return Vector3{x * other.x, y * other.y, z * other.z};
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
        }
        double std::min() const {
            const double t = abs(x) < abs(y) ? x : y;
            return abs(t) < abs(z) ? abs(t) : abs(z);
        }
        double aproject(const Vector3& u, const Vector3& v) const {
            Vector3 _cross = u.cross(v);
            if (_cross.len() == 0.0) return this->aproject(u);
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
