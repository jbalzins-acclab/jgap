#ifndef JGAP_VECTOR3_HPP
#define JGAP_VECTOR3_HPP

#include <cmath>
#include <string>


namespace jgap {

    /// Main structure for storing 3D vectors.
    struct Vector3 {
        constexpr Vector3() : x(0.0), y(0.0), z(0.0) {}
        constexpr Vector3(double x, double y, double z) : x(x), y(y), z(z) {}
        constexpr Vector3(const Vector3& other) = default;

        double x, y, z;

        constexpr Vector3 operator+(const Vector3& other) const {
            return Vector3{x + other.x, y + other.y, z + other.z};
        }

        constexpr Vector3& operator+=(const Vector3& other) {
            x += other.x;
            y += other.y;
            z += other.z;
            return *this;
        }

        constexpr Vector3 operator-() const {
            return Vector3{-x, -y, -z};
        }

        constexpr Vector3 operator-(const Vector3& other) const {
            return Vector3{x - other.x, y - other.y, z - other.z};
        }

        constexpr Vector3& operator-=(const Vector3& other) {
            x -= other.x;
            y -= other.y;
            z -= other.z;
            return *this;
        }

        constexpr Vector3 operator*(const double scalar) const { return Vector3{x * scalar, y * scalar, z * scalar}; }

        constexpr friend Vector3 operator*(const double scalar, const Vector3& vec) {
            return Vector3{vec.x * scalar, vec.y * scalar, vec.z * scalar};
        }

        constexpr Vector3 componentMul(const Vector3 other) const {
            return Vector3{x * other.x, y * other.y, z * other.z};
        }

        constexpr Vector3& operator*=(const double scalar) {
            x *= scalar;
            y *= scalar;
            z *= scalar;
            return *this;
        }

        constexpr Vector3 operator/(const double scalar) const { return Vector3{x / scalar, y / scalar, z / scalar}; }

        constexpr Vector3& operator/=(const double scalar) {
            x /= scalar;
            y /= scalar;
            z /= scalar;
            return *this;
        }

        constexpr double dot(const Vector3& other) const { return x * other.x + y * other.y + z * other.z; }

        constexpr Vector3 cross(const Vector3& other) const {
            return Vector3{y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x};
        }

        constexpr double square() const { return x * x + y * y + z * z; }

        double norm() const { return sqrt(x * x + y * y + z * z); }

        double project(const Vector3& other) const { return dot(other) / other.norm(); }

        /// @return This vector's component perpendicular to the other.
        double aproject(const Vector3& other) const { return sqrt(norm() * norm() - project(other) * project(other)); }

        /// @return This vector's component perpendicular to the plane in which u and v lie.
        double aproject(const Vector3& u, const Vector3& v) const {
            Vector3 _cross = u.cross(v);
            if (_cross.norm() == 0.0) [[unlikely]]
                return this->aproject(u);
            [[likely]] return abs(this->project(_cross));
        }

        Vector3 normalize() const { return *this * (1.0 / norm()); }

        double minComponent() const {
            const double t = abs(x) < abs(y) ? x : y;
            return abs(t) < abs(z) ? abs(t) : abs(z);
        }

        constexpr bool operator==(const Vector3& other) const { return x == other.x && y == other.y && z == other.z; }

        constexpr bool operator<(const Vector3& other) const {
            if (x < other.x) return true;
            if (x > other.x) return false;
            if (y < other.y) return true;
            if (y > other.y) return false;
            return z < other.z;
        }

        constexpr bool isPositive() const { return Vector3(0.0, 0.0, 0.0) < *this; }

        std::string toString() const { return std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z); }
    };
}

#endif
