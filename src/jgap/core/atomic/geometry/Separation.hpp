#ifndef JGAP_SEPARATION_HPP
#define JGAP_SEPARATION_HPP

#include "../../Vector3.hpp"
#include "jgap/core/Real.hpp"
#include "jgap/core/atomic/energy/Virials.hpp"

namespace jgap {

    /// @brief $|r_{ij}|$, $\hat{r}_{ij}$, and calculating @ref Virials of $|r_{ij}|$.
    ///
    /// Provides the best balance in terms of avoiding re-calculation and optimizing memory usage
    /// when storing positional information between a pair of atoms.
    ///
    /// @note To avoid any confusion down the line, unless performance critical simplification exists,
    /// Separation::virials() should be used as a starting point for the chain rule application
    /// when calculating Virials for e.g. @ref AtomicQuantity or @ref ManyBodyDescriptor.
    struct Separation {
        Real magnitude{};
        Vector3 direction{};

        Separation() = default;
        Separation(const Vector3 &r_from, const Vector3 &r_to) {
            const auto separation = r_to - r_from;
            magnitude = separation.norm();
            direction = separation / magnitude;
        }

        /// @brief $\grad_{\vec{r}_j} |r_{ij}|$
        Vector3 vec() const { return direction * magnitude; }

        /// @brief $\partial |r_{ij}| / \partial H_{ab}$ (@see Virials)
        Virials virials() const {
            return Virials{
                .xx = direction.x * direction.x,
                .xy = direction.x * direction.y,
                .xz = direction.x * direction.z,
                .yy = direction.y * direction.y,
                .yz = direction.y * direction.z,
                .zz = direction.z * direction.z
            } * -magnitude;
        }
    };
}

#endif