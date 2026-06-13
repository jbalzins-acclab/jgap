#include "Atoms.hpp"
#include <variant>
#include "io/log/CurrentLogger.hpp"
#include <sstream>

namespace jgap {
    Atoms::Atoms(const std::vector<Vector3>& pos,
                 const std::vector<Species>& spec,
                 const std::optional<Lattice>& lat,
                 std::array<bool, 3> pbc,
                 const AtomsPropertyNames& names)
        : XYZData(), main_property_names(names) {

        arrays[main_property_names.positions] = pos;
        arrays[main_property_names.species] = spec;

        if (lat) {
            properties[main_property_names.lattice] = *lat;
        }
        properties[main_property_names.pbc] = pbc;

        if (pos.size() != spec.size()) {
            JGAP_LOG_AND_THROW("Positions and species must have same size");
        }
        validateSizes();
        wrapPositions();
    }

    Atoms::Atoms(const XYZData& data, const AtomsPropertyNames& names)
        : XYZData(data), main_property_names(names) {

        bool has_pos = arrays.contains(main_property_names.positions)
            && std::holds_alternative<std::vector<Vector3>>(arrays.at(main_property_names.positions));
        bool has_species = arrays.contains(main_property_names.species)
            && std::holds_alternative<std::vector<Species>>(arrays.at(main_property_names.species));

        if (!has_pos || !has_species) {
            std::ostringstream oss;
            oss << "Atoms constructor from XYZData failed. Missing: ";
            if (!has_pos) oss << main_property_names.positions << " ";
            if (!has_species) oss << main_property_names.species << " ";
            oss << "\nAvailable arrays: ";
            for (auto const& [k, v] : arrays) oss << k << " ";
            JGAP_LOG_AND_THROW(oss.str());
        }

        validateSizes();
        wrapPositions();
    }

    std::optional<Lattice> Atoms::lookupLattice() const {
        if (properties.contains(main_property_names.lattice) && std::holds_alternative<Lattice>(properties.at(main_property_names.lattice))) {
            return std::get<Lattice>(properties.at(main_property_names.lattice));
        }
        return std::nullopt;
    }

    void Atoms::setLattice(const std::optional<Lattice>& lat) {
        if (lat) {
            properties[main_property_names.lattice] = *lat;
        } else {
            properties.erase(main_property_names.lattice);
        }
    }

    std::array<bool, 3> Atoms::lookupPbc() const {
        if (properties.contains(main_property_names.pbc)
            && std::holds_alternative<std::array<bool, 3>>(properties.at(main_property_names.pbc))) {
            return std::get<std::array<bool, 3>>(properties.at(main_property_names.pbc));
        }
        return {false, false, false};
    }

    void Atoms::setPbc(const std::array<bool, 3>& pbc) {
        properties[main_property_names.pbc] = pbc;
    }

    const std::vector<Vector3>& Atoms::lookupPositions() const {
        if (!arrays.contains(main_property_names.positions)
            || !std::holds_alternative<std::vector<Vector3>>(arrays.at(main_property_names.positions))) {
            JGAP_LOG_AND_THROW("Positions were deleted or altered incorrectly");
            }
        return std::get<std::vector<Vector3>>(arrays.at(main_property_names.positions));
    }

    std::vector<Vector3>& Atoms::lookupPositions() {
        if (!arrays.contains(main_property_names.positions)
            || !std::holds_alternative<std::vector<Vector3>>(arrays.at(main_property_names.positions))) {
            JGAP_LOG_AND_THROW("Positions were deleted or altered incorrectly");
        }
        return std::get<std::vector<Vector3>>(arrays[main_property_names.positions]);
    }

    const std::vector<Species>& Atoms::lookupSpecies() const {
        if (!arrays.contains(main_property_names.species)
            || !std::holds_alternative<std::vector<Species>>(arrays.at(main_property_names.species))) {
            JGAP_LOG_AND_THROW("Species were deleted or altered incorrectly");
            }
        return std::get<std::vector<Species>>(arrays.at(main_property_names.species));
    }

    std::vector<Species>& Atoms::lookupSpecies() {
        if (!arrays.contains(main_property_names.species)
            || !std::holds_alternative<std::vector<Species>>(arrays.at(main_property_names.species))) {
            JGAP_LOG_AND_THROW("Species were deleted or altered incorrectly");
            }
        return std::get<std::vector<Species>>(arrays[main_property_names.species]);
    }

    void Atoms::addAtom(const std::map<std::string, AtomValue>& atom_data) {
        if (!atom_data.contains(main_property_names.positions) || !atom_data.contains(main_property_names.species)) {
            JGAP_LOG_AND_THROW("addAtom requires positions and species arrays");
        }

        for (auto& [name, array] : arrays) {
            if (!atom_data.contains(name)) {
                JGAP_LOG_AND_THROW("Missing data for array: " + name + " when adding atom");
            }

            std::visit([&](auto&& arg) {
                using T = std::decay_t<decltype(arg)>;
                using ItemT = typename T::value_type;
                if (!std::holds_alternative<ItemT>(atom_data.at(name))) {
                    JGAP_LOG_AND_THROW("Type mismatch for array: " + name + " when adding atom");
                }
                arg.push_back(std::get<ItemT>(atom_data.at(name)));
            }, array);
        }
    }

    void Atoms::removeAtom(size_t index) {
        for (auto& [name, array] : arrays) {
            std::visit([index](auto&& arg) {
                arg.erase(arg.begin() + index);
            }, array);
        }
    }

    void Atoms::removeArray(const std::string& name) {
        if (name == main_property_names.positions || name == main_property_names.species) {
            JGAP_LOG_AND_THROW("Cannot remove core arrays: pos and species");
        }
        arrays.erase(name);
    }

    std::optional<Real> Atoms::getEnergy() const {
        if (properties.contains(main_property_names.energy)) {
            return std::get<Real>(properties.at(main_property_names.energy));
        }
        return std::nullopt;
    }

    void Atoms::lookupEnergy(Real e) {
        properties[main_property_names.energy] = e;
    }

    std::optional<Virials> Atoms::lookupVirials() const {
        if (properties.contains(main_property_names.virials)) {
            return std::get<Virials>(properties.at(main_property_names.virials));
        }
        return std::nullopt;
    }

    void Atoms::setVirials(const Virials& v) {
        properties[main_property_names.virials] = v;
    }

    std::optional<std::vector<Vector3>> Atoms::lookupForces() const {
        if (arrays.contains(main_property_names.forces)) {
            return std::get<std::vector<Vector3>>(arrays.at(main_property_names.forces));
        }
        return std::nullopt;
    }

    void Atoms::setForces(const std::vector<Vector3>& f) {
        if (f.size() != nAtoms()) {
            JGAP_LOG_AND_THROW("Forces size must match number of atoms");
        }
        arrays[main_property_names.forces] = f;
    }

    std::optional<std::string> Atoms::lookupConfigType() const {
        if (properties.contains(main_property_names.config_type) && std::holds_alternative<std::string>(properties.at(main_property_names.config_type))) {
            return std::get<std::string>(properties.at(main_property_names.config_type));
        }
        return std::nullopt;
    }

    void Atoms::setConfigType(const std::string& config_type) {
        properties[main_property_names.config_type] = config_type;
    }

    void Atoms::validateSizes() const {
        size_t n = lookupPositions().size();
        for (const auto& [name, array] : arrays) {
            size_t current_size = 0;
            std::visit([&current_size](auto&& arg) { current_size = arg.size(); }, array);
            if (n != current_size) {
                JGAP_LOG_AND_THROW("Size mismatch in Atoms: " + name + " has size " + std::to_string(current_size) + " but expected " + std::to_string(n));
            }
        }
    }

    void Atoms::wrapPositions() {
        auto pbc = lookupPbc();
        if (!pbc[0] && !pbc[1] && !pbc[2]) return;

        auto lat_opt = lookupLattice();
        if (!lat_opt) {
            JGAP_LOG_AND_THROW("PBC is true but no lattice is provided.");
        }

        const auto& lat = *lat_opt;
        Real V = lat.volume();

        if (std::abs(V) < 1e-12) {
            JGAP_LOG_AND_THROW("Lattice volume is too small or zero.");
        }

        // Calculate inverse matrix rows (reciprocal lattice vectors without 2pi)
        Vector3 r0 = lat.b.cross(lat.c) / V;
        Vector3 r1 = lat.c.cross(lat.a) / V;
        Vector3 r2 = lat.a.cross(lat.b) / V;

        for (auto& pos: lookupPositions()) {
            // Convert to fractional coordinates
            Real f0 = pos.dot(r0);
            Real f1 = pos.dot(r1);
            Real f2 = pos.dot(r2);

            // Wrap fractional coordinates to [0, 1) if PBC is enabled for that dimension
            if (pbc[0]) f0 -= std::floor(f0);
            if (pbc[1]) f1 -= std::floor(f1);
            if (pbc[2]) f2 -= std::floor(f2);

            // Convert back to Cartesian
            pos = lat.a * f0 + lat.b * f1 + lat.c * f2;
        }
    }
}