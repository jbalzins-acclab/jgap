#ifndef JGAP_SEPARATION_HPP
#define JGAP_SEPARATION_HPP

#include "Vector3.hpp"
#include "core/Real.hpp"
#include "core/atomic/energy/Virials.hpp"

namespace jgap {
    struct Separation {
        Real magnitude{};
        Vector3 direction{};
        Virials virials{};

        Separation() = default;
        Separation(const Separation&) = default;

        Separation(const Vector3 &r_from, const Vector3 &r_to) {
            const auto separation = r_to - r_from;
            magnitude = separation.norm();
            direction = separation.normalize();
            virials = {
                separation.x * direction.x,
                separation.x * direction.y,
                separation.x * direction.z,
                separation.y * direction.y,
                separation.y * direction.z,
                separation.z * direction.z
            };
            virials *= -1.0;
        }
        Vector3 vec() const { return direction * magnitude; }

        /*
        Virials virials() const {
            return Virials{
                direction.x * direction.x,
                direction.x * direction.y,
                direction.x * direction.z,
                direction.y * direction.y,
                direction.y * direction.z,
                direction.z * direction.z
            } * magnitude;
        }
        */
    };
}

#endif