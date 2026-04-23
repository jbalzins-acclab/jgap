#ifndef JGAP_VIRIALS_HPP
#define JGAP_VIRIALS_HPP

#include "core/Real.hpp"

namespace jgap {
    struct Virials {
        Real xx, xy, xz, yy, yz, zz;

        Virials operator+(const Virials& other) const {
            return {xx + other.xx, xy + other.xy, xz + other.xz,
                    yy + other.yy, yz + other.yz, zz + other.zz};
        }

        Virials operator-(const Virials& other) const {
            return {xx - other.xx, xy - other.xy, xz - other.xz,
                    yy - other.yy, yz - other.yz, zz - other.zz};
        }

        Virials& operator+=(const Virials& other) {
            xx += other.xx; xy += other.xy; xz += other.xz;
            yy += other.yy; yz += other.yz; zz += other.zz;
            return *this;
        }

        Virials& operator-=(const Virials& other) {
            xx -= other.xx; xy -= other.xy; xz -= other.xz;
            yy -= other.yy; yz -= other.yz; zz -= other.zz;
            return *this;
        }

        Virials operator*(Real scalar) const {
            return {xx * scalar, xy * scalar, xz * scalar,
                    yy * scalar, yz * scalar, zz * scalar};
        }

        Virials& operator*=(Real scalar) {
            xx *= scalar; xy *= scalar; xz *= scalar;
            yy *= scalar; yz *= scalar; zz *= scalar;
            return *this;
        }

        Virials& operator/=(Real scalar) {
            xx /= scalar; xy /= scalar; xz /= scalar;
            yy /= scalar; yz /= scalar; zz /= scalar;
            return *this;
        }
    };
}

#endif