#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "jgap/core/ValuePtr.hpp"
#include "jgap/core/Vector3.hpp"
#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/core/atomic/geometry/Lattice.hpp"
#include "jgap/core/atomic/species/Species.hpp"
#include "jgap/core/potentials/Potential.hpp"
#include "jgap/io/tabgap/TabGapIO.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"

namespace py = pybind11;
using namespace jgap;

class PyJGAPCalculator {
public:
    PyJGAPCalculator(const std::vector<std::string>& paths) {
        if (paths.empty()) {
            throw std::invalid_argument("Empty paths list");
        }

        bool is_tabgap = false;
        if (paths[0].length() >= 10 && paths[0].substr(paths[0].length() - 10) == ".tabgap.h5") {
            is_tabgap = true;
        }

        if (is_tabgap) {
            potential = TabGapIO::read(paths);
        } else {
            potential = SerializationRegistry<Potential>::deserialize(paths[0]);
        }
    }

    py::dict calculate(
        const std::vector<std::string>& chemical_symbols, py::array_t<double> positions, py::array_t<double> cell,
        py::array_t<bool> pbc
    ) {

        auto pos_buf = positions.request();
        auto cell_buf = cell.request();
        auto pbc_buf = pbc.request();

        size_t n_atoms = chemical_symbols.size();

        double* pos_ptr = static_cast<double*>(pos_buf.ptr);
        double* cell_ptr = static_cast<double*>(cell_buf.ptr);
        bool* pbc_ptr = static_cast<bool*>(pbc_buf.ptr);

        std::vector<Vector3> pos;
        pos.reserve(n_atoms);
        std::vector<Species> spec;
        spec.reserve(n_atoms);

        for (size_t i = 0; i < n_atoms; ++i) {
            pos.emplace_back(pos_ptr[3 * i], pos_ptr[3 * i + 1], pos_ptr[3 * i + 2]);
            spec.emplace_back(chemical_symbols[i]);
        }

        Lattice lat{
            Vector3(cell_ptr[0], cell_ptr[1], cell_ptr[2]), Vector3(cell_ptr[3], cell_ptr[4], cell_ptr[5]),
            Vector3(cell_ptr[6], cell_ptr[7], cell_ptr[8])
        };

        std::array<bool, 3> pbc_arr = {pbc_ptr[0], pbc_ptr[1], pbc_ptr[2]};

        Atoms atoms(pos, spec, lat, pbc_arr);

        AtomicQuantity result = potential->calculateEnergy(atoms);

        py::array_t<double> forces_array({n_atoms, (size_t) 3});
        auto forces_buf = forces_array.request();
        double* forces_ptr = static_cast<double*>(forces_buf.ptr);

        for (size_t i = 0; i < n_atoms; ++i) {
            forces_ptr[3 * i] = result.forces[i].x;
            forces_ptr[3 * i + 1] = result.forces[i].y;
            forces_ptr[3 * i + 2] = result.forces[i].z;
        }

        // Voigt order for ASE: xx, yy, zz, yz, xz, xy
        py::array_t<double> virial_array(6);
        auto virial_buf = virial_array.request();
        double* virial_ptr = static_cast<double*>(virial_buf.ptr);

        virial_ptr[0] = result.virials.xx;
        virial_ptr[1] = result.virials.yy;
        virial_ptr[2] = result.virials.zz;
        virial_ptr[3] = result.virials.yz;
        virial_ptr[4] = result.virials.xz;
        virial_ptr[5] = result.virials.xy;

        py::dict out;
        out["energy"] = result.value;
        out["forces"] = forces_array;
        out["virial"] = virial_array;
        return out;
    }

private:
    ValuePtr<Potential> potential;
};

PYBIND11_MODULE(jgap_ase, m) {
    py::class_<PyJGAPCalculator>(m, "PyJGAPCalculator")
        .def(py::init<const std::vector<std::string>&>())
        .def("calculate", &PyJGAPCalculator::calculate);
}
