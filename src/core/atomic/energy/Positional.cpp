#include "Positional.hpp"

#include <cassert>

namespace jgap {
    Positional Positional::operator+(const Positional &other) const {
        assert(forces.size() == other.forces.size() && "Number of forces doesn't match");

        Positional result;

        result.value = value + other.value;
        result.virials = virials + other.virials;

        result.forces = forces;
        for (size_t i = 0; i < forces.size(); i++) {
            result.forces[i] += other.forces[i];
        }

        return result;
    }

    Positional& Positional::operator+=(const Positional &other) {
        assert(forces.size() == other.forces.size() && "Number of forces doesn't match");

        value += other.value;
        virials += other.virials;

        for (size_t i = 0; i < forces.size(); i++) {
            forces[i] += other.forces[i];
        }

        return *this;
    }

    Positional Positional::operator*(double scalar) const {
        Positional result;

        result.value = value * scalar;
        result.virials = virials * scalar;

        result.forces = forces;
        for (size_t i = 0; i < forces.size(); i++) {
            result.forces[i] *= scalar;
        }

        return result;
    }

    Positional& Positional::operator*=(double scalar) {
        value *= scalar;
        virials *= scalar;
        for (auto& force : forces) {
            force *= scalar;
        }
        return *this;
    }

    // Chain rule - apply cutoff
    Positional Positional::operator*(const Positional &other) const {
        assert(forces.size() == other.forces.size() && "Number of forces doesn't match");

        Positional result;
        result.value = value * other.value;

        result.virials = virials * other.value + other.virials * value;

        result.forces.resize(forces.size());
        for (size_t i = 0; i < forces.size(); ++i) {
            result.forces[i] = forces[i] * other.value + other.forces[i] * value;
        }

        return result;

    }
}
