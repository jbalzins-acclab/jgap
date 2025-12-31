#include "PredictionData.hpp"

namespace jgap {
    Predictions Predictions::operator+(const Predictions &other) const {
        Predictions result;
        result.energy = energy + other.energy;
        result.virials = virials + other.virials;

        // Assuming forces vector sizes are compatible for addition
        // If one is empty and the other is not, the non-empty one is taken.
        // If both are empty, result.forces remains empty.
        // If both are non-empty, they are added element-wise.
        if (hasForces() && other.hasForces()) {
            result.forces_optional.resize(forces_optional.size());
            for (size_t i = 0; i < forces_optional.size(); ++i) {
                result.forces_optional[i] = forces_optional[i] + other.forces_optional[i];
            }
        } else if (!forces_optional.empty()) {
            result.forces_optional = forces_optional;
        } else if (!other.forces_optional.empty()) {
            result.forces_optional = other.forces_optional;
        }
        return result;
    }

    Predictions & Predictions::operator+=(const Predictions &other) {
        energy += other.energy;
        virials += other.virials;

        if (hasForces() && other.hasForces()) {
            for (size_t i = 0; i < forces_optional.size(); ++i) {
                forces_optional[i] += other.forces_optional[i];
            }
        } else if (!other.forces_optional.empty()) {
            // If current object has no forces but other does, copy other's forces
            forces_optional = other.forces_optional;
        }
        // If current object has forces and other does not, keep current forces.
        // If both are empty, they remain empty.
        return *this;
    }

    Predictions Predictions::operator*(double scalar) const {
        Predictions result;
        result.energy = energy * scalar;
        result.virials = virials * scalar;
        result.forces_optional.resize(forces_optional.size());
        for (size_t i = 0; i < forces_optional.size(); ++i) {
            result.forces_optional[i] = forces_optional[i] * scalar;
        }
        return result;
    }

    Predictions & Predictions::operator*=(double scalar) {
        energy *= scalar;
        virials *= scalar;
        for (auto& force : forces_optional) {
            force *= scalar;
        }
        return *this;
    }
}
