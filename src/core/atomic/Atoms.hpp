#ifndef JGAP_ATOMICBOX_HPP
#define JGAP_ATOMICBOX_HPP

#include <array>
#include <cmath>
#include <map>
#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "../Vector3.hpp"
#include "core/Real.hpp"
#include "energy/AtomicQuantity.hpp"
#include "geometry/Lattice.hpp"
#include "io/XYZData.hpp"
#include "species/Species.hpp"
#include "utils/Utils.hpp"

namespace jgap {

    /// Types underlying \ref XYZArrayType vector types.
    using PerAtomProperty = std::variant<int, Real, Vector3, std::string, Species>;

    /// Stores the positions and species of atoms in a structure,
    /// as well as information on periodicity, and may optionally store
    /// lattice vectors (will provide if the system is periodic)
    /// and energy/virials/forces,
    /// and the "config type" of the structure.
    ///
    /// Positions are guaranteed to be wrapped within the box (i.e. fractional coords \in [0; 1) )
    /// in periodic dimensions.
    ///
    /// \note Inspired by Atoms from ASE.
    // TODO: Update documentation to reflect decoupling from XYZData
    class Atoms {
    public:
        static std::vector<Atoms> readAtoms(const std::string& filename, const MainXYZPropertyNames& prop_names = {}) {
            return mapVector(XYZData::read(filename, prop_names), [](const XYZData& d) { return Atoms(d); });
        }

        Atoms(const std::vector<Vector3>& pos, const std::vector<Species>& spec,
              const std::optional<Lattice>& lat = std::nullopt, std::array<bool, 3> pbc = {false, false, false},
              MainXYZPropertyNames names = {});

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

        XYZData toXYZData() const;
        void write(const std::string& filename) const;
        void write(std::ostream& out_stream) const;
        std::string write() const;

        /// Set energy, forces and virials.
        /// @note essentially to simplify updating energies with the result of \ref Potential::calculateEnergy.
        void setEnergyAndDerivatives(const AtomicQuantity& energy_and_derivatives);
        void setEnergyAndDerivatives(Real energy);

        std::optional<Lattice> getLattice() const { return lattice; }
        void setLattice(const Lattice& lat) { lattice = lat; }
        void eraseLattice() { lattice = std::nullopt; }

        std::array<bool, 3> getPbc() const { return pbc; }
        void setPbc(const std::array<bool, 3>& p) { pbc = p; }

        const std::vector<Vector3>& getPositions() const { return positions; }
        std::vector<Vector3>& getPositions() { return positions; }

        const std::vector<Species>& getSpecies() const { return species; }
        std::vector<Species>& getSpecies() { return species; }

        size_t nAtoms() const { return positions.size(); }

        void addAtom(const std::map<std::string, PerAtomProperty>& atom_data);
        void removeAtom(size_t index);
        void removeArray(const std::string& name);

        std::optional<Real> getEnergy() const { return energy; }
        void setEnergy(Real e) { energy = e; }
        void eraseEnergy() { energy = std::nullopt; }

        std::optional<Virials> getVirials() const { return virials; }
        void setVirials(const Virials& v) { virials = v; }
        void eraseVirials() { virials = std::nullopt; }

        std::optional<std::vector<Vector3>> getForces() const { return forces; }
        void setForces(const std::vector<Vector3>& f) { forces = f; }
        void eraseForces() { forces = std::nullopt; }

        std::optional<std::string> getConfigType() const { return config_type; }
        void setConfigType(const std::string& type) { config_type = type; }
        void eraseConfigType() { config_type = std::nullopt; }

        const std::map<std::string, XYZInfoType>& getExtraInfo() const { return extra_info; }
        std::map<std::string, XYZInfoType>& getExtraInfoForEditing() { return extra_info; }

        const std::map<std::string, XYZArrayType>& getExtraArrays() const { return extra_arrays; }
        std::map<std::string, XYZArrayType>& getExtraArraysForEditing() { return extra_arrays; }

        const MainXYZPropertyNames& getMainPropertyNames() const { return main_property_names; }

        void wrapPositions();

    private:
        MainXYZPropertyNames main_property_names;
        std::vector<Vector3> positions;
        std::vector<Species> species;
        std::optional<std::vector<Vector3>> forces;
        std::optional<Lattice> lattice;
        std::array<bool, 3> pbc = {false, false, false};
        std::optional<Real> energy;
        std::optional<Virials> virials;
        std::optional<std::string> config_type;

        std::map<std::string, XYZInfoType> extra_info;
        std::map<std::string, XYZArrayType> extra_arrays;

        void validateXYZ(const XYZData& data);
        void validateSizes() const;
    };
}

#endif
