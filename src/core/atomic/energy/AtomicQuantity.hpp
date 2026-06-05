#ifndef JGAP_ATOMICQUANTITY_HPP
#define JGAP_ATOMICQUANTITY_HPP

#include "../geometry/Vector3.hpp"
#include <vector>

#include "Virials.hpp"

namespace jgap {
    struct AtomicQuantity {
        // itself
        Real value;
        // derivatives
        Virials virials;
        std::vector<Vector3> forces;

        AtomicQuantity() : value(Real{}), virials({}), forces({}) {}
        AtomicQuantity(size_t n_atoms) : value(Real{}), virials({}), forces(n_atoms, Vector3{}) {}

        bool empty() const {
            return value == 0.0 && forces.empty();
        }

        AtomicQuantity operator+(const AtomicQuantity& other) const;
        AtomicQuantity& operator+=(const AtomicQuantity& other);
        AtomicQuantity operator*(double scalar) const;
        AtomicQuantity& operator*=(double scalar);
    };
}

#endif
