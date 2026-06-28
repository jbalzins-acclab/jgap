#ifndef JGAP_CLUSTER_HPP
#define JGAP_CLUSTER_HPP

#include "../../Vector3.hpp"
#include "Separation.hpp"
#include <cassert>

#include "core/CalculationType.hpp"
#include "core/atomic/neighbours/NeighbourData.hpp"

namespace jgap {

    /// @brief Calculates flat-array index for a symmetric matrix.
    /// @note The first matrix index must be lower than the second index.
    /// To save performance, this is unchecked unless the DEBUG flag is on.
    static constexpr size_t flattenedIndex(size_t lower_index, size_t higher_index) {
#ifdef DEBUG
        assert(lower_index < higher_index);
#endif
        return higher_index * (higher_index - 1) / 2 + lower_index;
    }

    /// @brief Separation matrix for a cluster of atoms,
    /// as well as their indexes in Atoms from which they originate.
    /// \tparam NAtoms Number of atoms in a cluster.
    /// \tparam CalcType Indicates whether the \ref SeparationDerivatives should be stored as well.
    ///
    /// @note Only independent matrix components are stored in a flat array,
    /// indexed as in \ref flattenedIndex.
    template<size_t NAtoms, CalculationType CalcType = ValueOnly>
    requires(NAtoms > 1)
    struct Cluster;

    template<size_t NAtoms>
    struct Cluster<NAtoms, ValueOnly> {
        static constexpr size_t NSeparations = NAtoms * (NAtoms - 1) / 2;

        std::array<size_t, NAtoms> atom_indexes;
        std::array<Real, NSeparations> separation_magnitudes;

        Cluster() = default;
        Cluster(const Cluster&) = default;

        Cluster(size_t index0, const std::array<NeighbourData, NAtoms - 1>& atom_neigh) {

            atom_indexes[0] = index0;

            for (size_t i = 1; i < NAtoms; i++) {
                atom_indexes[i] = atom_neigh[i - 1].neighbour_index;
            }

            for (size_t j = 1; j < NAtoms; j++) {
                separation_magnitudes[flattenedIndex(0, j)] = atom_neigh[j - 1].separation.magnitude;
            }

            for (size_t i = 1; i < NAtoms; i++) {
                for (size_t j = i + 1; j < NAtoms; j++) {
                    separation_magnitudes[flattenedIndex(i, j)]
                        = (atom_neigh[i - 1].separation.vec() - atom_neigh[j - 1].separation.vec()).norm();
                }
            }
        }

        /// \see flattenedIndex
        Real separationBetween(const size_t lower_index, const size_t higher_index) const {
#ifdef DEBUG
            assert(lower_index < higher_index);
#endif
            return separation_magnitudes[flattenedIndex(lower_index, higher_index)];
        }

        /// \see flattenedIndex
        Real& separationBetween(const size_t lower_index, const size_t higher_index) {
#ifdef DEBUG
            assert(lower_index < higher_index);
#endif
            return separation_magnitudes[flattenedIndex(lower_index, higher_index)];
        }
    };

    template<size_t NAtoms>
    struct Cluster<NAtoms, WithGradients> : Cluster<NAtoms> {
        using Base = Cluster<NAtoms>;
        using Base::NSeparations;
        using Base::atom_indexes;
        using Base::separation_magnitudes;

        std::array<SeparationDerivatives, NSeparations> derivatives;

        Cluster() = default;
        Cluster(const Cluster&) = default;

        Cluster(size_t index0, const std::array<NeighbourData, NAtoms - 1>& atom_neigh) {

            atom_indexes[0] = index0;

            for (size_t i = 1; i < NAtoms; i++) {
                atom_indexes[i] = atom_neigh[i - 1].neighbour_index;
            }

            for (size_t j = 1; j < NAtoms; j++) {
                separation_magnitudes[flattenedIndex(0, j)] = atom_neigh[j - 1].separation.magnitude;
                derivatives[flattenedIndex(0, j)] = atom_neigh[j - 1].separation.derivatives;
            }

            for (size_t i = 1; i < NAtoms; i++) {
                for (size_t j = i + 1; j < NAtoms; j++) {

                    Vector3 relative_pos_i = atom_neigh[i - 1].separation.vec();
                    Vector3 relative_pos_j = atom_neigh[j - 1].separation.vec();

                    auto separation = Separation(relative_pos_i, relative_pos_j);

                    separation_magnitudes[flattenedIndex(i, j)] = separation.magnitude;
                    derivatives[flattenedIndex(i, j)] = separation.derivatives;
                }
            }
        }

        /// \see flattenedIndex
        const SeparationDerivatives& derivativesBetween(const size_t lower_index, const size_t higher_index) const {
#ifdef DEBUG
            assert(lower_index < higher_index);
#endif
            return derivatives[flattenedIndex(lower_index, higher_index)];
        }

        /// \see flattenedIndex
        SeparationDerivatives& derivativesBetween(const size_t lower_index, const size_t higher_index) {
#ifdef DEBUG
            assert(lower_index < higher_index);
#endif
            return derivatives[flattenedIndex(lower_index, higher_index)];
        }
    };
}

#endif