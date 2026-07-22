#ifndef JGAP_ATOMICQUANTITIES_HPP
#define JGAP_ATOMICQUANTITIES_HPP

#include <cassert>
#include <vector>
#include "AtomicQuantity.hpp"
#include "core/Matrix.hpp"

namespace jgap {

    /// \brief An array of \ref AtomicQuantity, but with a tighter memory layout.
    /// \note Avoids having a vector of force vectors when dealing with per-sparse point kernel values.
    class AtomicQuantities {
    public:
        AtomicQuantities(size_t n_sparse, size_t n_atoms) :
            n_sparse(n_sparse), n_atoms(n_atoms), energy_data(n_sparse, Real{}), virial_data(n_sparse, Virials{}),
            force_table(n_sparse * n_atoms, Vector3{}) {}

        Real &energy(size_t sparse_idx) { return energy_data[sparse_idx]; }
        const Real &energy(size_t sparse_idx) const { return energy_data[sparse_idx]; }

        Virials &virials(size_t sparse_idx) { return virial_data[sparse_idx]; }
        const Virials &virials(size_t sparse_idx) const { return virial_data[sparse_idx]; }

        Vector3 &force(size_t sparse_idx, size_t atom_idx) { return force_table[sparse_idx * n_atoms + atom_idx]; }
        const Vector3 &force(size_t sparse_idx, size_t atom_idx) const {
            return force_table[sparse_idx * n_atoms + atom_idx];
        }

        AtomicQuantity reduce(const std::vector<Real> &coefficients) const {
            assert(coefficients.size() == n_sparse);

            AtomicQuantity total(n_atoms);
            for (size_t i = 0; i < n_sparse; ++i) {
                total.value += energy_data[i] * coefficients[i];
                total.virials += virial_data[i] * coefficients[i];
                for (size_t j = 0; j < n_atoms; j++) {
                    total.forces[j] += force(i, j) * coefficients[i];
                }
            }
            return total;
        }

        AtomicQuantities& operator*=(Real scalar) {
            for (auto& e : energy_data) e *= scalar;
            for (auto& v : virial_data) v *= scalar;
            for (auto& f : force_table) f *= scalar;
            return *this;
        }

    private:
        size_t n_sparse;
        size_t n_atoms;

        std::vector<Real> energy_data;
        std::vector<Virials> virial_data;
        std::vector<Vector3> force_table;
    };
}

#endif
