#ifndef JGAP_POTENTIALPREDICTION_HPP
#define JGAP_POTENTIALPREDICTION_HPP

#include "../Vector3.hpp"
#include <vector>
#include <optional>

namespace jgap {

    struct Virials {
        double xx, xy, xz, yy, yz, zz;

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

        Virials operator*(double scalar) const {
            return {xx * scalar, xy * scalar, xz * scalar,
                    yy * scalar, yz * scalar, zz * scalar};
        }

        Virials& operator*=(double scalar) {
            xx *= scalar; xy *= scalar; xz *= scalar;
            yy *= scalar; yz *= scalar; zz *= scalar;
            return *this;
        }

        Virials& operator/=(double scalar) {
            xx /= scalar; xy /= scalar; xz /= scalar;
            yy /= scalar; yz /= scalar; zz /= scalar;
            return *this;
        }
    };

    struct Predictions {
        double energy;
        Virials virials;
        std::vector<Vector3> forces_optional;

        Predictions operator+(const Predictions& other) const;
        Predictions& operator+=(const Predictions& other);
        Predictions operator*(double scalar) const;
        Predictions& operator*=(double scalar);

        bool hasForces() const { return !forces_optional.empty(); }
    };
}

#endif
