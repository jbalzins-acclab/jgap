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

namespace jgap {

    using AtomValue = std::variant<int, Real, Vector3, std::string, Species>;

    struct AtomsPropertyNames {
        std::string positions = "pos";
        std::string species = "species";
        std::string forces = "force";
        std::string virials = "virial";
        std::string energy = "energy";
        std::string lattice = "Lattice";
        std::string pbc = "pbc";
        std::string config_type = "config_type";
    };

    class Atoms : public XYZData {
    public:
        Atoms(const std::vector<Vector3> &pos,
              const std::vector<Species> &spec,
              const std::optional<Lattice> &lat = std::nullopt,
              std::array<bool, 3> pbc = {false, false, false},
              const AtomsPropertyNames &names = {});

        explicit Atoms(const XYZData& data, const AtomsPropertyNames& names = {});

        Atoms(const Atoms& other) = default;
        Atoms& operator=(const Atoms& other) = default;

        Atoms& operator<<(const AtomicQuantity& quantity);

        std::optional<Lattice> lookupLattice() const;
        void setLattice(const std::optional<Lattice>& lat);

        std::array<bool, 3> lookupPbc() const;
        void setPbc(const std::array<bool, 3>& pbc);

        const std::vector<Vector3>& lookupPositions() const;
        std::vector<Vector3>& lookupPositions();

        const std::vector<Species>& lookupSpecies() const;
        std::vector<Species>& lookupSpecies();

        size_t nAtoms() const { return lookupPositions().size(); }

        void addAtom(const std::map<std::string, AtomValue>& atom_data);
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

        const AtomsPropertyNames& getMainPropertyNames() const { return main_property_names; }

    private:
        AtomsPropertyNames main_property_names;

        void validateSizes() const;
        void wrapPositions();
    };
}

#endif