#include "Atoms.hpp"
#include <utility>
#include <variant>
#include "../io/log/CurrentLogger.hpp"
#include <sstream>

namespace jgap {
    Atoms::Atoms(const std::vector<Vector3>& pos,
                 const std::vector<Species>& spec,
                 const std::optional<Lattice>& lat,
                 std::array<bool, 3> pbc_val,
                 MainXYZPropertyNames names)
        : main_property_names(std::move(names)), positions(pos), species(spec), lattice(lat), pbc(pbc_val) {

        if (pos.size() != spec.size()) {
            JGAP_LOG_AND_THROW("Positions and species must have same size");
        }
        validateSizes();
        wrapPositions();
    }

    Atoms::Atoms(const XYZData& data) : main_property_names(data.getMainPropertyNames()) {
        validateXYZ(data);
    }

    Atoms::Atoms(XYZData &&data) : main_property_names(data.getMainPropertyNames()) {
        validateXYZ(data);
    }

    XYZData Atoms::toXYZData() const {
        XYZData data(main_property_names);
        auto& info = data.getPropertiesForEditing();
        auto& arrs = data.getArraysForEditing();

        // Copy extra properties first
        info = extra_info;
        arrs = extra_arrays;

        // Overwrite core properties
        arrs[main_property_names.positions] = positions;
        arrs[main_property_names.species] = species;

        if (lattice) {
            info[main_property_names.lattice] = *lattice;
        }
        info[main_property_names.pbc] = pbc;

        if (energy) {
            info[main_property_names.energy] = *energy;
        }
        if (virials) {
            info[main_property_names.virials] = *virials;
        }
        if (forces) {
            arrs[main_property_names.forces] = *forces;
        }
        if (config_type) {
            info[main_property_names.config_type] = *config_type;
        }

        return data;
    }

    void Atoms::write(const std::string &filename) const {
        toXYZData().write(filename);
    }

    void Atoms::write(std::ostream &out_stream) const {
        toXYZData().write(out_stream);
    }

    std::string Atoms::write() const {
        return toXYZData().write();
    }

    void Atoms::setEnergyAndDerivatives(const AtomicQuantity& energy_and_derivatives) {
        setEnergy(energy_and_derivatives.value);
        setVirials(energy_and_derivatives.virials);
        setForces(energy_and_derivatives.forces);
    }

    void Atoms::addAtom(const std::map<std::string, PerAtomProperty>& atom_data) {
        if (!atom_data.contains(main_property_names.positions) || !atom_data.contains(main_property_names.species)) {
            JGAP_LOG_AND_THROW("addAtom requires positions and species arrays");
        }

        positions.push_back(std::get<Vector3>(atom_data.at(main_property_names.positions)));
        species.push_back(std::get<Species>(atom_data.at(main_property_names.species)));

        if (forces) {
            if (!atom_data.contains(main_property_names.forces)) {
                JGAP_LOG_AND_THROW("Missing data for forces when adding atom");
            }
            forces->push_back(std::get<Vector3>(atom_data.at(main_property_names.forces)));
        }

        for (auto& [name, array] : extra_arrays) {
            if (!atom_data.contains(name)) {
                JGAP_LOG_AND_THROW("Missing data for array: {} when adding atom", name);
            }

            std::visit([&](auto&& arg) {
                using T = std::decay_t<decltype(arg)>;
                using ItemT = typename T::value_type;
                if (!std::holds_alternative<ItemT>(atom_data.at(name))) {
                    JGAP_LOG_AND_THROW("Type mismatch for array: {} when adding atom", name);
                }
                arg.push_back(std::get<ItemT>(atom_data.at(name)));
            }, array);
        }
    }

    void Atoms::removeAtom(size_t index) {
        positions.erase(positions.begin() + index);
        species.erase(species.begin() + index);
        
        if (forces) {
            forces->erase(forces->begin() + index);
        }

        for (auto& [name, array] : extra_arrays) {
            std::visit([index](auto&& arg) {
                arg.erase(arg.begin() + index);
            }, array);
        }
    }

    void Atoms::removeArray(const std::string& name) {
        if (name == main_property_names.positions || name == main_property_names.species) {
            JGAP_LOG_AND_THROW("Cannot remove core arrays: pos and species");
        }
        if (name == main_property_names.forces) {
            forces = std::nullopt;
            return;
        }
        extra_arrays.erase(name);
    }

    void Atoms::validateXYZ(const XYZData& data) {
        const auto& arrays = data.getArrays();
        const auto& info = data.getProperties();

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

        positions = std::get<std::vector<Vector3>>(arrays.at(main_property_names.positions));
        species = std::get<std::vector<Species>>(arrays.at(main_property_names.species));

        if (info.contains(main_property_names.lattice) && std::holds_alternative<Lattice>(info.at(main_property_names.lattice))) {
            lattice = std::get<Lattice>(info.at(main_property_names.lattice));
        }

        if (info.contains(main_property_names.pbc) && std::holds_alternative<std::array<bool, 3>>(info.at(main_property_names.pbc))) {
            pbc = std::get<std::array<bool, 3>>(info.at(main_property_names.pbc));
        }

        if (info.contains(main_property_names.energy) && std::holds_alternative<Real>(info.at(main_property_names.energy))) {
            energy = std::get<Real>(info.at(main_property_names.energy));
        }

        if (info.contains(main_property_names.virials) && std::holds_alternative<Virials>(info.at(main_property_names.virials))) {
            virials = std::get<Virials>(info.at(main_property_names.virials));
        }

        if (arrays.contains(main_property_names.forces) && std::holds_alternative<std::vector<Vector3>>(arrays.at(main_property_names.forces))) {
            forces = std::get<std::vector<Vector3>>(arrays.at(main_property_names.forces));
        }

        if (info.contains(main_property_names.config_type) && std::holds_alternative<std::string>(info.at(main_property_names.config_type))) {
            config_type = std::get<std::string>(info.at(main_property_names.config_type));
        }

        // Copy extra arrays
        for (const auto& [name, array] : arrays) {
            if (name != main_property_names.positions &&
                name != main_property_names.species &&
                name != main_property_names.forces) {
                extra_arrays[name] = array;
            }
        }

        // Copy extra info
        for (const auto& [name, val] : info) {
            if (name != main_property_names.lattice &&
                name != main_property_names.pbc &&
                name != main_property_names.energy &&
                name != main_property_names.virials &&
                name != main_property_names.config_type) {
                extra_info[name] = val;
            }
        }

        validateSizes();
        wrapPositions();
    }

    void Atoms::validateSizes() const {
        size_t n = positions.size();
        if (species.size() != n) {
            JGAP_LOG_AND_THROW("Size mismatch in Atoms: species has size {} but expected {}", species.size(), n);
        }
        if (forces && forces->size() != n) {
            JGAP_LOG_AND_THROW("Size mismatch in Atoms: forces has size {} but expected {}", forces->size(), n);
        }

        for (const auto& [name, array] : extra_arrays) {
            size_t current_size = 0;
            std::visit([&current_size](auto&& arg) { current_size = arg.size(); }, array);
            if (n != current_size) {
                JGAP_LOG_AND_THROW("Size mismatch in Atoms: {} has size {} but expected {}", name, current_size, n);
            }
        }
    }

    void Atoms::wrapPositions() {
        if (!pbc[0] && !pbc[1] && !pbc[2]) return;

        if (!lattice) {
            JGAP_LOG_AND_THROW("PBC is true but no lattice is provided.");
        }

        const auto& lat = *lattice;
        Real V = lat.volume();

        if (std::abs(V) < 1e-12) {
            JGAP_LOG_AND_THROW("Lattice volume is too small or zero.");
        }

        // Rows of matrix that transforms coordinates s.t. lattice vectors are orthonormal.
        Vector3 r0 = lat.b.cross(lat.c) / V; // => .a = 1
        Vector3 r1 = lat.c.cross(lat.a) / V; // => .b = 1
        Vector3 r2 = lat.a.cross(lat.b) / V; // => .c = 1

        for (auto& pos: positions) {
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