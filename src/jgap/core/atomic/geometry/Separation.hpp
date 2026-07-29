#ifndef JGAP_SEPARATION_HPP
#define JGAP_SEPARATION_HPP

#include "../../Vector3.hpp"
#include "jgap/core/Real.hpp"
#include "jgap/core/atomic/energy/Virials.hpp"

namespace jgap {

    /// @brief $|r_{ij}|$ and derivatives.
    ///
    /// Calculating virial stress tensor of $|r_{ij}|$ here ensures
    /// no convention related confusion further down the line.
    ///
    /// @note Direction is practically always needed as without it one can encounter problems
    /// related to PBC in 2-body symmetric iteration.
    ///
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