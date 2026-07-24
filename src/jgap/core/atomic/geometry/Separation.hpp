#ifndef JGAP_SEPARATION_HPP
#define JGAP_SEPARATION_HPP

#include "../../Vector3.hpp"
#include "jgap/core/Real.hpp"
#include "jgap/core/atomic/energy/Virials.hpp"

namespace jgap {

    /// @brief $\grad_{\vec{r}_j} |r_{ij}|$ and $\partial |r_{ij}| / \partial H_ab$ (\see Virials).
    struct SeparationDerivatives {
        Vector3 direction{};
        Virials virials{};

        SeparationDerivatives() = default;
        SeparationDerivatives(const SeparationDerivatives&) = default;
        SeparationDerivatives(SeparationDerivatives&&) = default;
        SeparationDerivatives& operator=(const SeparationDerivatives&) = default;
        SeparationDerivatives& operator=(SeparationDerivatives&&) = default;
    };

    /// @brief $|r_{ij}|$ and derivatives(\see SeparationDerivatives).
    ///
    /// Calculating virial stress tensor of $|r_{ij}|$ here ensures
    /// no convention related confusion further down the line.
    ///
    /// @note Storing the original separation vector and re-calculating magnitude/direction each time is slower
    /// than one SIMD friendly vec() = (direction * magnitude).
    /// Storing both vectors adds strain on memory, which is also slower.
    /// @note Direction is practically always needed as without it one can encounter problems
    /// related to PBC in 2-body symmetric iteration.
    /// @note Virial calculation on demand showed no performance enhancement compared with pre-calculating;
    /// however, the test was not quite systematic, so one might re-try if interested.
    ///
    struct Separation {
        Real magnitude{};
        SeparationDerivatives derivatives{};

        Separation() = default;
        Separation(const Vector3 &r_from, const Vector3 &r_to) {
            const auto separation = r_to - r_from;
            magnitude = separation.norm();
            derivatives.direction = separation / magnitude;
            derivatives.virials = {
                separation.x * derivatives.direction.x,
                separation.x * derivatives.direction.y,
                separation.x * derivatives.direction.z,
                separation.y * derivatives.direction.y,
                separation.y * derivatives.direction.z,
                separation.z * derivatives.direction.z
            };
            derivatives.virials *= -1.0;
        }
        Vector3 vec() const { return derivatives.direction * magnitude; }
    };
}

#endif