#ifndef JGAP_ATOMICQUANTITY_HPP
#define JGAP_ATOMICQUANTITY_HPP

#include "../../Vector3.hpp"
#include <vector>

#include "Virials.hpp"

namespace jgap {
    /// Quantity (like energy, or kernel value) based on atomic coordinates,
    /// and it's derivatives wrt to them and a linear transformation on them.
    ///
    /// @note Essentially equivalent to a one-dimensional \ref ManyBodyDescriptor,
    /// except fixed Dim here allows avoiding bypass the use of \ref std::array
    /// which is redundant here, for no-template calculations.
    ///
    struct AtomicQuantity {
        Real value;

        Virials virials;
        std::vector<Vector3> forces;

        explicit AtomicQuantity(size_t n_atoms) : value(Real{}), virials({}), forces(n_atoms, Vector3{}) {}

        AtomicQuantity operator+(const AtomicQuantity& other) const;
        AtomicQuantity& operator+=(const AtomicQuantity& other);
        AtomicQuantity operator*(Real scalar) const;
        AtomicQuantity& operator*=(Real scalar);
    };
}

#endif
