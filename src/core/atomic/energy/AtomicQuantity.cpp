#include "AtomicQuantity.hpp"

#include <cassert>

namespace jgap {
    AtomicQuantity AtomicQuantity::operator+(const AtomicQuantity &other) const {
        if (empty()) return other;
        if (other.empty()) return *this;

        assert(forces.size() == other.forces.size() && "Number of forces doesn't match");

        AtomicQuantity result(forces.size());

        result.value = value + other.value;
        result.virials = virials + other.virials;

        result.forces = forces;
        for (size_t i = 0; i < forces.size(); i++) {
            result.forces[i] += other.forces[i];
        }

        return result;
    }

    AtomicQuantity& AtomicQuantity::operator+=(const AtomicQuantity &other) {
        if (other.empty()) return *this;
        if (empty()) {
            *this = other;
            return *this;
        }

        assert(forces.size() == other.forces.size() && "Number of forces doesn't match");

        value += other.value;
        virials += other.virials;

        for (size_t i = 0; i < forces.size(); i++) {
            forces[i] += other.forces[i];
        }

        return *this;
    }

    AtomicQuantity AtomicQuantity::operator*(double scalar) const {
        AtomicQuantity result(forces.size());

        result.value = value * scalar;
        result.virials = virials * scalar;

        result.forces = forces;
        for (size_t i = 0; i < forces.size(); i++) {
            result.forces[i] *= scalar;
        }

        return result;
    }

    AtomicQuantity& AtomicQuantity::operator*=(double scalar) {
        value *= scalar;
        virials *= scalar;
        for (auto& force : forces) {
            force *= scalar;
        }
        return *this;
    }
}
