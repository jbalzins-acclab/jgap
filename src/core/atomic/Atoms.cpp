#include "Atoms.hpp"
#include <variant>
#include "io/log/CurrentLogger.hpp"
#include <sstream>

namespace jgap {

    Atoms::Atoms(const std::vector<Vector3>& pos, const std::vector<Species>& spec, const std::optional<Lattice>& lat, std::array<bool, 3> pbc)
        : positions(std::get<std::vector<Vector3>>(arrays["pos"] = pos)),
          species(std::get<std::vector<Species>>(arrays["species"] = spec)) {
        if (lat) {
            properties["Lattice"] = *lat;
        }
        properties["pbc"] = pbc;

        if (pos.size() != spec.size()) {
            JGAP_LOG_AND_THROW("Positions and species must have same size");
        }
        validateSizes();
    }

    Atoms::Atoms(const XYZData& data)
        : XYZData(data),
          positions(std::get<std::vector<Vector3>>(arrays.at("pos"))),
          species(std::get<std::vector<Species>>(arrays.at("species"))) {

        bool has_pos = arrays.contains("pos") && std::holds_alternative<std::vector<Vector3>>(arrays.at("pos"));
        bool has_species = arrays.contains("species") && std::holds_alternative<std::vector<Species>>(arrays.at("species"));

        if (!has_pos || !has_species) {
            std::ostringstream oss;
            oss << "Atoms constructor from XYZData failed. Missing: ";
            if (!has_pos) oss << "pos ";
            if (!has_species) oss << "species ";
            oss << "\nAvailable arrays: ";
            for (auto const& [k, v] : arrays) oss << k << " ";
            JGAP_LOG_AND_THROW(oss.str());
        }
        validateSizes();
    }

    std::optional<Lattice> Atoms::getLattice() const {
        if (properties.contains("Lattice") && std::holds_alternative<Lattice>(properties.at("Lattice"))) {
            return std::get<Lattice>(properties.at("Lattice"));
        }
        return std::nullopt;
    }

    void Atoms::setLattice(const std::optional<Lattice>& lat) {
        if (lat) {
            properties["Lattice"] = *lat;
        } else {
            properties.erase("Lattice");
        }
    }

    std::array<bool, 3> Atoms::getPbc() const {
        if (properties.contains("pbc") && std::holds_alternative<std::array<bool, 3>>(properties.at("pbc"))) {
            return std::get<std::array<bool, 3>>(properties.at("pbc"));
        }
        return {false, false, false};
    }

    void Atoms::setPbc(const std::array<bool, 3>& pbc) {
        properties["pbc"] = pbc;
    }

    void Atoms::addAtom(const std::map<std::string, AtomValue>& atom_data) {
        if (!atom_data.contains("pos") || !atom_data.contains("species")) {
            JGAP_LOG_AND_THROW("addAtom requires 'pos' and 'species'");
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
        if (index >= positions.size()) {
            JGAP_LOG_AND_THROW("Index out of bounds in removeAtom");
        }

        for (auto& [name, array] : arrays) {
            std::visit([index](auto&& arg) {
                arg.erase(arg.begin() + index);
            }, array);
        }
    }

    void Atoms::removeArray(const std::string& name) {
        if (name == "pos" || name == "species") {
            JGAP_LOG_AND_THROW("Cannot remove core arrays: pos and species");
        }
        arrays.erase(name);
    }

    std::optional<Real> Atoms::getEnergy() const {
        if (properties.contains("energy")) {
            return std::get<Real>(properties.at("energy"));
        }
        return std::nullopt;
    }

    void Atoms::setEnergy(Real e) {
        properties["energy"] = e;
    }

    std::optional<Virials> Atoms::getVirials() const {
        if (properties.contains("virial")) {
            return std::get<Virials>(properties.at("virial"));
        }
        return std::nullopt;
    }

    void Atoms::setVirials(const Virials& v) {
        properties["virial"] = v;
    }

    std::optional<std::vector<Vector3>> Atoms::getForces() const {
        if (arrays.contains("force")) {
            return std::get<std::vector<Vector3>>(arrays.at("force"));
        }
        return std::nullopt;
    }

    void Atoms::setForces(const std::vector<Vector3>& f) {
        if (f.size() != nAtoms()) {
            JGAP_LOG_AND_THROW("Forces size must match number of atoms");
        }
        arrays["force"] = f;
    }

    std::optional<std::string> Atoms::getConfigType() const {
        if (properties.contains("config_type") && std::holds_alternative<std::string>(properties.at("config_type"))) {
            return std::get<std::string>(properties.at("config_type"));
        }
        return std::nullopt;
    }

    void Atoms::setConfigType(const std::string& config_type) {
        properties["config_type"] = config_type;
    }

    void Atoms::validateSizes() const {
        size_t n = positions.size();
        for (const auto& [name, array] : arrays) {
            size_t current_size = 0;
            std::visit([&current_size](auto&& arg) { current_size = arg.size(); }, array);
            if (n != current_size) {
                JGAP_LOG_AND_THROW("Size mismatch in Atoms: " + name + " has size " + std::to_string(current_size) + " but expected " + std::to_string(n));
            }
        }
    }
}
