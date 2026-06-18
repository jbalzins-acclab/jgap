#ifndef JGAP_SEPARATION_HPP
#define JGAP_SEPARATION_HPP

#include "Vector3.hpp"
#include "core/Real.hpp"
#include "core/atomic/energy/Virials.hpp"

namespace jgap {

    struct SeparationDerivatives {
        Vector3 direction{};
        Virials virials{};

        SeparationDerivatives() = default;
        SeparationDerivatives(const SeparationDerivatives&) = default;
    };

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