#ifndef JGAP_ATOMICBOX_HPP
#define JGAP_ATOMICBOX_HPP

#include <optional>
#include <map>
#include <vector>
#include <string>
#include <variant>
#include <array>
#include <cmath>

#include "../Vector3.hpp"
#include "species/Species.hpp"
#include "geometry/Lattice.hpp"
#include "core/Real.hpp"
#include "io/XYZData.hpp"
#include "energy/AtomicQuantity.hpp"
#include "utils/Utils.hpp"

namespace jgap {

    /// Types underlying \ref XYZArrayType vector types.
    using PerAtomPropery = std::variant<int, Real, Vector3, std::string, Species>;

    /// Stores the positions and species of atoms in a structure,
    /// as well as information on periodicity, and may optionally store
    /// lattice vectors (will provide if the system is periodic)
    /// and energy/virials/forces,
    /// and the "config type" of the structure.
    ///
    /// Positions are guaranteed to be wrapped within the box (i.e. fractional coords \in [0; 1) )
    /// in periodic dimensions.
    ///
    /// \note Derived from \ref XYZData to avoid losing any relevant data from the input files.
    /// To avoid duplication and errors that would come with it,
    /// all quantities are found by searching them in the \ref XYZData's maps,
    /// looking up by property names defined in \ref MainXYZPropertyNames upon construction.
    /// This comes at the cost of repeated access being inefficient,
    /// so the getters use the prefix "lookup" instead of "get"
    /// to indicate that finding property requires a map lookup.
    ///
    /// \note Inspired by Atoms from ASE.
    class Atoms : public XYZData {
    public:

        static std::vector<Atoms> readAtoms(const std::string& filename, const MainXYZPropertyNames& prop_names = {}) {
            return mapVector(read(filename, prop_names), [&prop_names](const XYZData& d) { return Atoms(d); });
        }

        Atoms(const std::vector<Vector3> &pos,
              const std::vector<Species> &spec,
              const std::optional<Lattice> &lat = std::nullopt,
              std::array<bool, 3> pbc = {false, false, false},
              const MainXYZPropertyNames &names = {});

        /// Construct by coping \ref XYZData.
        /// \throws \ref std::runtime_error if \ref XYZData doesn't contain
        /// position/species arrays disguised by array name defined in \ref MainXYZPropertyNames,
        /// or if PBC info in XYZData is inconsistent.
        explicit Atoms(const XYZData& data);

        /// Construct by coping \ref XYZData.
        /// \throws \ref std::runtime_error if \ref XYZData doesn't contain
        /// position/species arrays disguised by array name defined in \ref MainXYZPropertyNames,
        /// or if PBC info in XYZData is inconsistent.
        explicit Atoms(XYZData&& data);

        Atoms(const Atoms& other) = default;
        Atoms(Atoms&& other) = default;
        Atoms& operator=(const Atoms& other) = default;
        Atoms& operator=(Atoms&& other) = default;

        /// Set energy, forces and virials.
        /// @note essentially to simplify updating energies with the result of \ref Potential::calculateEnergy.
        Atoms& operator<<(const AtomicQuantity& energy_and_derivatives);

        std::optional<Lattice> lookupLattice() const;
        void setLattice(const std::optional<Lattice>& lat);

        std::array<bool, 3> lookupPbc() const;
        void setPbc(const std::array<bool, 3>& pbc);

        const std::vector<Vector3>& lookupPositions() const;
        std::vector<Vector3>& lookupPositions();

        const std::vector<Species>& lookupSpecies() const;
        std::vector<Species>& lookupSpecies();

        size_t nAtoms() const { return lookupPositions().size(); }

        void addAtom(const std::map<std::string, PerAtomPropery>& atom_data);
        void removeAtom(size_t index);
        void removeArray(const std::string& name);

        std::optional<Real> lookupEnergy() const;
        void setEnergy(Real e);

        std::optional<Virials> lookupVirials() const;
        void setVirials(const Virials& v);

        std::optional<std::vector<Vector3>> lookupForces() const;
        void setForces(const std::vector<Vector3>& f);

        std::optional<std::string> lookupConfigType() const;
        void setConfigType(const std::string& config_type);

        void validateXYZ();
        void validateSizes() const;
        void wrapPositions();
    };
}

#endif