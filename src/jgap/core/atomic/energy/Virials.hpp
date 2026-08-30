#ifndef JGAP_VIRIALS_HPP
#define JGAP_VIRIALS_HPP

#include "jgap/core/Real.hpp"
#include "jgap/core/Vector3.hpp"

namespace jgap {

    /// Independent components of a virial stress tensor of positional quantity, following convention: <br>
    /// $V_{ij}(f) = -1 * \partial f / \partial H_{ij}$ |(evaluated at) H = 3D Identity matrix,
    /// where H is a linear transformation applied to the system's coordinates,
    /// such that $\vec{r}_i = H\vec{r}^{(0)}_i$.
    ///
    /// @note No volume normalization in this convention (i.e., without 1/V factor).
    /// @warning To avoid all kinds of errors,
    /// it is strongly advised to avoid calculating it directly unless performance-critical.
    /// Instead, put coordinate pairs (r_i, r_j) into @ref Separation, which will calculate the virials of |r_ij|,
    /// and then apply the chain rule so that V(q) = \partial q / \partial |r_ij| * V(|r_ij|).
    ///
    struct Virials {
        Real xx{}, xy{}, xz{}, yy{}, yz{}, zz{};

        /// @brief Constructs the symmetric virial stress tensor from the outer product (dyadic) $\vec{r} \otimes
        /// \vec{f}$.
        static constexpr Virials dyadic(const Vector3& r, const Vector3& f) {
            return { r.x * f.x, r.x * f.y, r.x * f.z, r.y * f.y, r.y * f.z, r.z * f.z };
        }

        Virials operator+(const Virials& other) const {
            return { xx + other.xx, xy + other.xy, xz + other.xz, yy + other.yy, yz + other.yz, zz + other.zz };
        }

        Virials operator-(const Virials& other) const {
            return { xx - other.xx, xy - other.xy, xz - other.xz, yy - other.yy, yz - other.yz, zz - other.zz };
        }

        Virials& operator+=(const Virials& other) {
            xx += other.xx;
            xy += other.xy;
            xz += other.xz;
            yy += other.yy;
            yz += other.yz;
            zz += other.zz;
            return *this;
        }

        Virials& operator-=(const Virials& other) {
            xx -= other.xx;
            xy -= other.xy;
            xz -= other.xz;
            yy -= other.yy;
            yz -= other.yz;
            zz -= other.zz;
            return *this;
        }

        Virials operator*(Real scalar) const {
            return { xx * scalar, xy * scalar, xz * scalar, yy * scalar, yz * scalar, zz * scalar };
        }

        friend Virials operator*(Real scalar, const Virials& v) { return v * scalar; }

        Virials& operator*=(Real scalar) {
            xx *= scalar;
            xy *= scalar;
            xz *= scalar;
            yy *= scalar;
            yz *= scalar;
            zz *= scalar;
            return *this;
        }

        Virials& operator/=(Real scalar) {
            Real inv = 1.0 / scalar;
            xx *= inv;
            xy *= inv;
            xz *= inv;
            yy *= inv;
            yz *= inv;
            zz *= inv;
            return *this;
        }
    };
}

#endif
