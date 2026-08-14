#ifndef JGAP_ATOMICQUANTITY_HPP
#define JGAP_ATOMICQUANTITY_HPP

#include "../../Vector3.hpp"
#include <vector>

#include "Virials.hpp"

namespace jgap {

    /// @brief F(Atoms), $\grad_{\vec{r}_i}F$ for all i, and its @ref Virials.
    /// @note Essentially equivalent to a one-dimensional @ref ManyBodyDescriptor<Dim = 1>,
    /// but is meant to store energy or scalar multiple of it (e.g. kernel value).
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
