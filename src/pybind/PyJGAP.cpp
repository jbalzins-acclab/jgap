#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "jgap/core/Vector3.hpp"
#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/core/atomic/energy/AtomicQuantity.hpp"
#include "jgap/core/atomic/energy/Virials.hpp"
#include "jgap/core/atomic/geometry/Lattice.hpp"
#include "jgap/core/atomic/species/Species.hpp"
#include "jgap/core/fit/gap/regularization/RegularizationRules.hpp"
#include "jgap/core/fit/gap/regularization/SimpleRegularizationRules.hpp"
#include "jgap/core/potentials/Cutoffs.hpp"
#include "jgap/core/potentials/Potential.hpp"
#include "jgap/core/potentials/gap/GapPotential.hpp"
#include "jgap/core/potentials/tabgap/TabGapPotential.hpp"
#include "jgap/core/tabulation/TabulationParams.hpp"
#include "jgap/core/transform/nbody/2b/eam/EamPairFunction.hpp"
#include "jgap/core/transform/nbody/2b/eam/FSGenPairFunction.hpp"
#include "jgap/io/PotentialLoader.hpp"
#include "jgap/io/tabgap/TabGapIO.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"
#include "jgap/utils/gap/StandardGapFit.hpp"
#include "jgap/utils/gap/StandardGapParams.hpp"

namespace py = pybind11;
using namespace jgap;

namespace {

    py::array_t<double> positionsToNumpy(const std::vector<Vector3>& pos) {
        size_t n = pos.size();
        py::array_t<double> arr({n, (size_t) 3});
        auto buf = arr.request();
        double* ptr = static_cast<double*>(buf.ptr);
        for (size_t i = 0; i < n; ++i) {
            ptr[3 * i + 0] = pos[i].x;
            ptr[3 * i + 1] = pos[i].y;
            ptr[3 * i + 2] = pos[i].z;
        }
        return arr;
    }

    std::vector<Vector3> numpyToPositions(py::array_t<double> arr) {
        auto buf = arr.request();
        if (buf.ndim == 2) {
            if (buf.shape[1] != 3) {
                throw std::invalid_argument("Expected 2D array with shape (N, 3)");
            }
            size_t n = buf.shape[0];
            const double* ptr = static_cast<const double*>(buf.ptr);
            std::vector<Vector3> pos;
            pos.reserve(n);
            for (size_t i = 0; i < n; ++i) {
                pos.emplace_back(ptr[3 * i + 0], ptr[3 * i + 1], ptr[3 * i + 2]);
            }
            return pos;
        } else if (buf.ndim == 1) {
            if (buf.shape[0] % 3 != 0) {
                throw std::invalid_argument("Expected 1D array with length multiple of 3");
            }
            size_t n = buf.shape[0] / 3;
            const double* ptr = static_cast<const double*>(buf.ptr);
            std::vector<Vector3> pos;
            pos.reserve(n);
            for (size_t i = 0; i < n; ++i) {
                pos.emplace_back(ptr[3 * i + 0], ptr[3 * i + 1], ptr[3 * i + 2]);
            }
            return pos;
        }
        throw std::invalid_argument("Expected 1D or 2D array for positions");
    }

    py::array_t<double> latticeToNumpy(const Lattice& lat) {
        py::array_t<double> arr({(size_t) 3, (size_t) 3});
        auto buf = arr.request();
        double* ptr = static_cast<double*>(buf.ptr);
        ptr[0] = lat.a.x; ptr[1] = lat.a.y; ptr[2] = lat.a.z;
        ptr[3] = lat.b.x; ptr[4] = lat.b.y; ptr[5] = lat.b.z;
        ptr[6] = lat.c.x; ptr[7] = lat.c.y; ptr[8] = lat.c.z;
        return arr;
    }

    Lattice numpyToLattice(py::array_t<double> arr) {
        auto buf = arr.request();
        const double* ptr = static_cast<const double*>(buf.ptr);
        if (buf.ndim == 2 && buf.shape[0] == 3 && buf.shape[1] == 3) {
            return Lattice{
                Vector3(ptr[0], ptr[1], ptr[2]),
                Vector3(ptr[3], ptr[4], ptr[5]),
                Vector3(ptr[6], ptr[7], ptr[8])
            };
        } else if (buf.ndim == 1 && buf.shape[0] == 9) {
            return Lattice{
                Vector3(ptr[0], ptr[1], ptr[2]),
                Vector3(ptr[3], ptr[4], ptr[5]),
                Vector3(ptr[6], ptr[7], ptr[8])
            };
        }
        throw std::invalid_argument("Expected (3, 3) or (9,) array for Lattice");
    }

    py::array_t<double> virialsToVoigt(const Virials& v) {
        py::array_t<double> arr(6);
        auto buf = arr.request();
        double* ptr = static_cast<double*>(buf.ptr);
        ptr[0] = v.xx;
        ptr[1] = v.yy;
        ptr[2] = v.zz;
        ptr[3] = v.yz;
        ptr[4] = v.xz;
        ptr[5] = v.xy;
        return arr;
    }

    py::array_t<double> virialsToMatrix(const Virials& v) {
        py::array_t<double> arr({(size_t) 3, (size_t) 3});
        auto buf = arr.request();
        double* ptr = static_cast<double*>(buf.ptr);
        ptr[0] = v.xx; ptr[1] = v.xy; ptr[2] = v.xz;
        ptr[3] = v.xy; ptr[4] = v.yy; ptr[5] = v.yz;
        ptr[6] = v.xz; ptr[7] = v.yz; ptr[8] = v.zz;
        return arr;
    }

    Virials voigtToVirials(py::array_t<double> arr) {
        auto buf = arr.request();
        const double* ptr = static_cast<const double*>(buf.ptr);
        if (buf.ndim == 1 && buf.shape[0] == 6) {
            return Virials{ptr[0], ptr[5], ptr[4], ptr[1], ptr[3], ptr[2]};
        } else if (buf.ndim == 2 && buf.shape[0] == 3 && buf.shape[1] == 3) {
            return Virials{ptr[0], ptr[1], ptr[2], ptr[4], ptr[5], ptr[8]};
        }
        throw std::invalid_argument("Expected (6,) Voigt array or (3, 3) matrix for Virials");
    }

    // Fast ASE calculation wrapper
    class PyJGAPCalculator {
    public:
        explicit PyJGAPCalculator(const std::vector<std::string>& paths) {
            potential = loadPotential(paths);
        }

        py::dict calculate(
            const std::vector<std::string>& chemical_symbols,
            py::array_t<double> positions,
            py::array_t<double> cell,
            py::array_t<bool> pbc
        ) {
            auto pos_buf = positions.request();
            auto cell_buf = cell.request();
            auto pbc_buf = pbc.request();

            size_t n_atoms = chemical_symbols.size();

            const double* pos_ptr = static_cast<const double*>(pos_buf.ptr);
            const double* cell_ptr = static_cast<const double*>(cell_buf.ptr);
            const bool* pbc_ptr = static_cast<const bool*>(pbc_buf.ptr);

            std::vector<Vector3> pos;
            pos.reserve(n_atoms);
            std::vector<Species> spec;
            spec.reserve(n_atoms);

            for (size_t i = 0; i < n_atoms; ++i) {
                pos.emplace_back(pos_ptr[3 * i], pos_ptr[3 * i + 1], pos_ptr[3 * i + 2]);
                spec.emplace_back(chemical_symbols[i]);
            }

            Lattice lat{
                Vector3(cell_ptr[0], cell_ptr[1], cell_ptr[2]),
                Vector3(cell_ptr[3], cell_ptr[4], cell_ptr[5]),
                Vector3(cell_ptr[6], cell_ptr[7], cell_ptr[8])
            };

            std::array<bool, 3> pbc_arr = {pbc_ptr[0], pbc_ptr[1], pbc_ptr[2]};

            Atoms atoms(pos, spec, lat, pbc_arr);

            AtomicQuantity result = potential->calculateEnergy(atoms);

            py::array_t<double> forces_array({n_atoms, (size_t) 3});
            auto forces_buf = forces_array.request();
            double* forces_ptr = static_cast<double*>(forces_buf.ptr);

            for (size_t i = 0; i < n_atoms; ++i) {
                forces_ptr[3 * i + 0] = result.forces[i].x;
                forces_ptr[3 * i + 1] = result.forces[i].y;
                forces_ptr[3 * i + 2] = result.forces[i].z;
            }

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

} // anonymous namespace

PYBIND11_MODULE(_jgap, m) {
    m.doc() = "jgap Python C++ bindings";

    // =========================================================================
    // Vector3
    // =========================================================================
    py::class_<Vector3>(m, "Vector3")
        .def(py::init<>())
        .def(py::init<Real, Real, Real>(), py::arg("x"), py::arg("y"), py::arg("z"))
        .def_readwrite("x", &Vector3::x)
        .def_readwrite("y", &Vector3::y)
        .def_readwrite("z", &Vector3::z)
        .def("norm", &Vector3::norm)
        .def("dot", &Vector3::dot)
        .def("cross", &Vector3::cross)
        .def("to_numpy", [](const Vector3& v) {
            py::array_t<double> arr(3);
            auto buf = arr.request();
            double* ptr = static_cast<double*>(buf.ptr);
            ptr[0] = v.x; ptr[1] = v.y; ptr[2] = v.z;
            return arr;
        })
        .def_static("from_numpy", [](py::array_t<double> arr) {
            auto buf = arr.request();
            if (buf.size != 3) throw std::invalid_argument("Expected array of size 3");
            const double* ptr = static_cast<const double*>(buf.ptr);
            return Vector3(ptr[0], ptr[1], ptr[2]);
        })
        .def("__add__", [](const Vector3& a, const Vector3& b) { return a + b; })
        .def("__sub__", [](const Vector3& a, const Vector3& b) { return a - b; })
        .def("__mul__", [](const Vector3& a, Real s) { return a * s; })
        .def("__rmul__", [](const Vector3& a, Real s) { return a * s; })
        .def("__truediv__", [](const Vector3& a, Real s) { return a / s; })
        .def("__repr__", [](const Vector3& v) {
            std::ostringstream oss;
            oss << "Vector3(" << v.x << ", " << v.y << ", " << v.z << ")";
            return oss.str();
        });

    // =========================================================================
    // Species
    // =========================================================================
    py::class_<Species>(m, "Species")
        .def(py::init<const std::string&>(), py::arg("symbol"))
        .def_static("from_atomic_number", &Species::fromAtomicNumber, py::arg("z"))
        .def_property_readonly("symbol", &Species::symbol)
        .def_property_readonly("id", &Species::getId)
        .def_property_readonly("atomic_number", &Species::atomicNumber)
        .def_property_readonly("mass", &Species::mass)
        .def("__repr__", [](const Species& s) {
            return "<Species '" + s.symbol() + "'>";
        })
        .def("__str__", &Species::symbol)
        .def("__eq__", &Species::operator==)
        .def("__hash__", [](const Species& s) { return std::hash<uint16_t>{}(s.getId()); });

    // =========================================================================
    // Lattice
    // =========================================================================
    py::class_<Lattice>(m, "Lattice")
        .def(py::init<>())
        .def(py::init<Vector3, Vector3, Vector3>(), py::arg("a"), py::arg("b"), py::arg("c"))
        .def(py::init(&numpyToLattice), py::arg("matrix"))
        .def_readwrite("a", &Lattice::a)
        .def_readwrite("b", &Lattice::b)
        .def_readwrite("c", &Lattice::c)
        .def("volume", &Lattice::volume)
        .def("to_numpy", &latticeToNumpy)
        .def("__repr__", [](const Lattice& lat) {
            std::ostringstream oss;
            oss << "Lattice(a=" << lat.a.x << "," << lat.a.y << "," << lat.a.z
                << ", b=" << lat.b.x << "," << lat.b.y << "," << lat.b.z
                << ", c=" << lat.c.x << "," << lat.c.y << "," << lat.c.z << ")";
            return oss.str();
        });

    // =========================================================================
    // Virials
    // =========================================================================
    py::class_<Virials>(m, "Virials")
        .def(py::init<>())
        .def(py::init<Real, Real, Real, Real, Real, Real>(),
             py::arg("xx"), py::arg("xy"), py::arg("xz"),
             py::arg("yy"), py::arg("yz"), py::arg("zz"))
        .def(py::init(&voigtToVirials), py::arg("voigt_or_matrix"))
        .def_readwrite("xx", &Virials::xx)
        .def_readwrite("xy", &Virials::xy)
        .def_readwrite("xz", &Virials::xz)
        .def_readwrite("yy", &Virials::yy)
        .def_readwrite("yz", &Virials::yz)
        .def_readwrite("zz", &Virials::zz)
        .def("to_voigt", &virialsToVoigt)
        .def("to_matrix", &virialsToMatrix)
        .def("__repr__", [](const Virials& v) {
            std::ostringstream oss;
            oss << "Virials(xx=" << v.xx << ", yy=" << v.yy << ", zz=" << v.zz
                << ", yz=" << v.yz << ", xz=" << v.xz << ", xy=" << v.xy << ")";
            return oss.str();
        });

    // =========================================================================
    // AtomicQuantity
    // =========================================================================
    py::class_<AtomicQuantity>(m, "AtomicQuantity")
        .def_readonly("value", &AtomicQuantity::value)
        .def_readonly("virials", &AtomicQuantity::virials)
        .def_property_readonly("energy", [](const AtomicQuantity& q) { return q.value; })
        .def_property_readonly("forces", [](const AtomicQuantity& q) {
            return positionsToNumpy(q.forces);
        })
        .def("__repr__", [](const AtomicQuantity& q) {
            std::ostringstream oss;
            oss << "<AtomicQuantity energy=" << q.value << ", n_forces=" << q.forces.size() << ">";
            return oss.str();
        });

    // =========================================================================
    // Atoms
    // =========================================================================
    py::class_<Atoms>(m, "Atoms")
        .def(py::init([](
            py::array_t<double> pos,
            const std::vector<std::string>& symbols,
            std::optional<Lattice> lat,
            std::optional<std::array<bool, 3>> pbc
        ) {
            auto pos_vec = numpyToPositions(pos);
            std::vector<Species> spec_vec;
            spec_vec.reserve(symbols.size());
            for (const auto& s : symbols) {
                spec_vec.emplace_back(s);
            }
            std::array<bool, 3> pbc_arr = pbc.value_or(std::array<bool, 3>{false, false, false});
            return Atoms(pos_vec, spec_vec, lat, pbc_arr);
        }),
        py::arg("positions"),
        py::arg("symbols"),
        py::arg("lattice") = std::nullopt,
        py::arg("pbc") = std::nullopt)

        .def(py::init([](
            const std::vector<Vector3>& pos,
            const std::vector<Species>& spec,
            std::optional<Lattice> lat,
            std::optional<std::array<bool, 3>> pbc
        ) {
            std::array<bool, 3> pbc_arr = pbc.value_or(std::array<bool, 3>{false, false, false});
            return Atoms(pos, spec, lat, pbc_arr);
        }),
        py::arg("positions"),
        py::arg("species"),
        py::arg("lattice") = std::nullopt,
        py::arg("pbc") = std::nullopt)

        .def_property("positions",
            [](const Atoms& a) { return positionsToNumpy(a.getPositions()); },
            [](Atoms& a, py::array_t<double> pos) { a.getPositions() = numpyToPositions(pos); })

        .def_property("species",
            [](const Atoms& a) { return a.getSpecies(); },
            [](Atoms& a, const std::vector<Species>& s) { a.getSpecies() = s; })

        .def_property("symbols",
            [](const Atoms& a) {
                std::vector<std::string> syms;
                syms.reserve(a.nAtoms());
                for (const auto& s : a.getSpecies()) syms.push_back(s.symbol());
                return syms;
            },
            [](Atoms& a, const std::vector<std::string>& syms) {
                std::vector<Species> spec;
                spec.reserve(syms.size());
                for (const auto& s : syms) spec.emplace_back(s);
                a.getSpecies() = spec;
            })

        .def_property("lattice",
            [](const Atoms& a) { return a.getLattice(); },
            [](Atoms& a, std::optional<Lattice> lat) {
                if (lat.has_value()) a.setLattice(*lat);
                else a.eraseLattice();
            })

        .def_property("pbc", &Atoms::getPbc, &Atoms::setPbc)

        .def_property("energy", &Atoms::getEnergy, [](Atoms& a, std::optional<Real> e) {
            if (e.has_value()) a.setEnergy(*e);
            else a.eraseEnergy();
        })

        .def_property("forces",
            [](const Atoms& a) -> std::optional<py::array_t<double>> {
                auto f = a.getForces();
                if (!f.has_value()) return std::nullopt;
                return positionsToNumpy(*f);
            },
            [](Atoms& a, std::optional<py::array_t<double>> f) {
                if (f.has_value()) a.setForces(numpyToPositions(*f));
                else a.eraseForces();
            })

        .def_property("virials", &Atoms::getVirials, [](Atoms& a, std::optional<Virials> v) {
            if (v.has_value()) a.setVirials(*v);
            else a.eraseVirials();
        })

        .def_property("config_type", &Atoms::getConfigType, [](Atoms& a, std::optional<std::string> ct) {
            if (ct.has_value()) a.setConfigType(*ct);
            else a.eraseConfigType();
        })

        .def("n_atoms", &Atoms::nAtoms)
        .def("__len__", &Atoms::nAtoms)
        .def("wrap_positions", &Atoms::wrapPositions)
        .def("write", py::overload_cast<const std::string&>(&Atoms::write, py::const_))

        .def_static("read_atoms", [](const std::string& filename) {
            return Atoms::readAtoms(filename);
        }, py::arg("filename"))

        .def("__repr__", [](const Atoms& a) {
            std::ostringstream oss;
            oss << "<Atoms n_atoms=" << a.nAtoms() << ", pbc=["
                << a.getPbc()[0] << ", " << a.getPbc()[1] << ", " << a.getPbc()[2] << "]>";
            return oss.str();
        });

    // =========================================================================
    // Cutoffs
    // =========================================================================
    py::class_<Cutoffs>(m, "Cutoffs")
        .def(py::init<>())
        .def("max_overall", &Cutoffs::maxOverall)
        .def("for_dim", &Cutoffs::forDim, py::arg("dim"))
        .def_readonly("per_cluster_size", &Cutoffs::per_cluster_size)
        .def("__repr__", [](const Cutoffs& c) {
            std::ostringstream oss;
            oss << "<Cutoffs max=" << c.maxOverall() << ">";
            return oss.str();
        });

    // =========================================================================
    // TabulationParams
    // =========================================================================
    py::class_<TabulationParams>(m, "TabulationParams")
        .def(py::init([](
            std::optional<Cutoffs> max_cutoffs,
            Real r_min_3b,
            Real max_eam_density,
            size_t n_grid_2b,
            std::optional<std::array<size_t, 3>> n_grid_3b
        ) {
            TabulationParams p;
            if (max_cutoffs.has_value()) p.max_cutoffs = *max_cutoffs;
            p.r_min_3b = r_min_3b;
            p.max_eam_density = max_eam_density;
            p.n_grid_2b = n_grid_2b;
            if (n_grid_3b.has_value()) p.n_grid_3b = *n_grid_3b;
            return p;
        }),
        py::arg("max_cutoffs") = std::nullopt,
        py::arg("r_min_3b") = 0.1,
        py::arg("max_eam_density") = 12.0,
        py::arg("n_grid_2b") = 5000,
        py::arg("n_grid_3b") = std::nullopt)
        .def_readwrite("max_cutoffs", &TabulationParams::max_cutoffs)
        .def_readwrite("r_min_3b", &TabulationParams::r_min_3b)
        .def_readwrite("max_eam_density", &TabulationParams::max_eam_density)
        .def_readwrite("n_grid_2b", &TabulationParams::n_grid_2b)
        .def_readwrite("n_grid_3b", &TabulationParams::n_grid_3b)
        .def("__repr__", [](const TabulationParams& p) {
            std::ostringstream oss;
            oss << "<TabulationParams r_min_3b=" << p.r_min_3b
                << ", max_eam_density=" << p.max_eam_density
                << ", n_grid_2b=" << p.n_grid_2b
                << ", n_grid_3b=[" << p.n_grid_3b[0] << "," << p.n_grid_3b[1] << "," << p.n_grid_3b[2] << "]>";
            return oss.str();
        });

    // =========================================================================
    // Potential
    // =========================================================================
    py::class_<Potential, std::shared_ptr<Potential>> pot_class(m, "Potential");

    // =========================================================================
    // GapPotential
    // =========================================================================
    py::class_<GapPotential, Potential, std::shared_ptr<GapPotential>>(m, "GapPotential")
        .def(py::init<>())
        .def("num_components", [](const GapPotential& p) { return p.getComponents().size(); })
        .def("__repr__", [](const GapPotential& p) {
            std::ostringstream oss;
            oss << "<GapPotential num_components=" << p.getComponents().size() << ">";
            return oss.str();
        });

    // =========================================================================
    // TabGapPotential
    // =========================================================================
    py::class_<TabGapPotential, Potential, std::shared_ptr<TabGapPotential>>(m, "TabGapPotential")
        .def(py::init<>())
        .def("__repr__", [](const TabGapPotential&) {
            return "<TabGapPotential>";
        });

    // Potential method bindings
    pot_class
        .def("calculate_energy", &Potential::calculateEnergy,
             py::arg("atoms"),
             py::call_guard<py::gil_scoped_release>(),
             "Calculate energy, forces, and virials for the given structure")
        .def("get_cutoffs", &Potential::getCutoffs)
        .def("tabulate", [](
            const Potential& pot,
            std::optional<TabulationParams> maybe_params,
            std::optional<Real> r_min_3b,
            std::optional<Real> max_eam_density,
            std::optional<size_t> n_grid_2b,
            std::optional<std::array<size_t, 3>> n_grid_3b
        ) -> std::shared_ptr<TabGapPotential> {
            TabulationParams params = maybe_params.value_or(TabulationParams{pot.getCutoffs()});
            if (r_min_3b.has_value()) params.r_min_3b = *r_min_3b;
            if (max_eam_density.has_value()) params.max_eam_density = *max_eam_density;
            if (n_grid_2b.has_value()) params.n_grid_2b = *n_grid_2b;
            if (n_grid_3b.has_value()) params.n_grid_3b = *n_grid_3b;

            TabulationData tab_data = pot.tabulate(params);
            return std::make_shared<TabGapPotential>(tab_data);
        },
        py::arg("params") = std::nullopt,
        py::arg("r_min_3b") = std::nullopt,
        py::arg("max_eam_density") = std::nullopt,
        py::arg("n_grid_2b") = std::nullopt,
        py::arg("n_grid_3b") = std::nullopt,
        py::call_guard<py::gil_scoped_release>(),
        "Tabulate potential into a TabGapPotential")
        .def("save", [](const Potential& pot, const std::string& filename) {
            SerializationRegistry<Potential>::serialize(
                ValuePtr<Potential>(std::unique_ptr<Potential>(pot.clone())), filename
            );
        }, py::arg("filename"), "Save potential to serialized HDF5 file")
        .def_static("load", [](const std::string& path) -> std::shared_ptr<Potential> {
            auto val = loadPotential(path);
            return std::shared_ptr<Potential>(val.release());
        }, py::arg("path"), "Load potential from file path or base stem");

    // =========================================================================
    // TabGapIO
    // =========================================================================
    py::class_<TabGapIO>(m, "TabGapIO")
        .def_static("write", [](const TabGapPotential& pot, const std::string& prefix) {
            return TabGapIO::write(pot, prefix);
        }, py::arg("potential"), py::arg("prefix"),
           "Write TabGapPotential to .tabgap.h5 and .eam.fs files")
        .def_static("read", [](const std::vector<std::string>& files) {
            return TabGapIO::read(files);
        }, py::arg("files"), "Read TabGapPotential from list of files")
        .def_static("read", [](const std::string& file) {
            return TabGapIO::read(std::vector<std::string>{file});
        }, py::arg("file"), "Read TabGapPotential from a single file");

    m.def("save_tabgap", [](const TabGapPotential& pot, const std::string& prefix) {
        return TabGapIO::write(pot, prefix);
    }, py::arg("potential"), py::arg("prefix"), "Save a TabGapPotential to disk");

    // =========================================================================
    // EAM Enums and Pair Functions
    // =========================================================================
    py::enum_<EamMode>(m, "EamMode")
        .value("FSsym", EamMode::FSsym)
        .value("FSgen", EamMode::FSgen)
        .value("EAM", EamMode::EAM)
        .value("Blind", EamMode::Blind)
        .export_values();

    py::class_<EamPairFunction, std::shared_ptr<EamPairFunction>>(m, "EamPairFunction");

    py::class_<FSGenPairFunction, EamPairFunction, std::shared_ptr<FSGenPairFunction>>(m, "FSGenPairFunction")
        .def(py::init<Real, Real, Real>(),
             py::arg("cutoff") = 4.5,
             py::arg("degree") = 3.0,
             py::arg("prefactor") = 1.0)
        .def_property_readonly("degree", &FSGenPairFunction::getDegree)
        .def("__repr__", [](const FSGenPairFunction& pf) {
            std::ostringstream oss;
            oss << "<FSGenPairFunction degree=" << pf.getDegree() << ">";
            return oss.str();
        });

    // =========================================================================
    // Regularization Rules
    // =========================================================================
    py::class_<RegularizationRules, std::shared_ptr<RegularizationRules>>(m, "RegularizationRules");

    py::class_<SimpleRegularizationRules, RegularizationRules, std::shared_ptr<SimpleRegularizationRules>>(m, "SimpleRegularizationRules")
        .def(py::init<Real, Real, Real, Real, Real, Real>(),
             py::arg("energy_sigma_per_atom") = 0.001,
             py::arg("force_component_sigma") = 0.05,
             py::arg("virials_iso_sigma_per_atom") = 0.1,
             py::arg("virials_aniso_sigmas_per_atom") = 0.02,
             py::arg("liquid_multiplier") = 5.0,
             py::arg("short_range_multiplier") = 5.0)
        .def("__repr__", [](const SimpleRegularizationRules&) {
            return "<SimpleRegularizationRules>";
        });

    // =========================================================================
    // StandardGapParams
    // =========================================================================
    py::class_<utils::StandardGapParams>(m, "StandardGapParams")
        .def(py::init([](
            size_t seed,
            std::optional<std::string> screened_coulomb_dataset_file,
            Real cutoff2,
            Real cutoff2_width,
            size_t n_sparse2,
            EamMode eam_mode,
            std::shared_ptr<EamPairFunction> eam_pf,
            size_t eam_n_sparse,
            Real eam_min_density,
            Real cutoff3,
            Real cutoff3_width,
            size_t n_sparse3,
            std::shared_ptr<RegularizationRules> regularization_rules
        ) {
            utils::StandardGapParams p;
            p.seed = seed;
            p.screened_coulomb_dataset_file = screened_coulomb_dataset_file;
            p.cutoff2 = cutoff2;
            p.cutoff2_width = cutoff2_width;
            p.n_sparse2 = n_sparse2;
            p.eam_mode = eam_mode;
            if (eam_pf) {
                p.eam_pf = ValuePtr<EamPairFunction>(std::unique_ptr<EamPairFunction>(eam_pf->clone()));
            }
            p.eam_n_sparse = eam_n_sparse;
            p.eam_min_density = eam_min_density;
            p.cutoff3 = cutoff3;
            p.cutoff3_width = cutoff3_width;
            p.n_sparse3 = n_sparse3;
            if (regularization_rules) {
                p.regularization_rules = ValuePtr<RegularizationRules>(std::unique_ptr<RegularizationRules>(regularization_rules->clone()));
            }
            return p;
        }),
        py::arg("seed") = 42,
        py::arg("screened_coulomb_dataset_file") = std::nullopt,
        py::arg("cutoff2") = 4.5,
        py::arg("cutoff2_width") = 1.0,
        py::arg("n_sparse2") = 20,
        py::arg("eam_mode") = EamMode::Blind,
        py::arg("eam_pf") = nullptr,
        py::arg("eam_n_sparse") = 20,
        py::arg("eam_min_density") = 0.05,
        py::arg("cutoff3") = 3.7,
        py::arg("cutoff3_width") = 0.6,
        py::arg("n_sparse3") = 500,
        py::arg("regularization_rules") = nullptr)
        .def_readwrite("seed", &utils::StandardGapParams::seed)
        .def_readwrite("screened_coulomb_dataset_file", &utils::StandardGapParams::screened_coulomb_dataset_file)
        .def_readwrite("cutoff2", &utils::StandardGapParams::cutoff2)
        .def_readwrite("cutoff2_width", &utils::StandardGapParams::cutoff2_width)
        .def_readwrite("n_sparse2", &utils::StandardGapParams::n_sparse2)
        .def_readwrite("eam_mode", &utils::StandardGapParams::eam_mode)
        .def_readwrite("eam_n_sparse", &utils::StandardGapParams::eam_n_sparse)
        .def_readwrite("eam_min_density", &utils::StandardGapParams::eam_min_density)
        .def_readwrite("cutoff3", &utils::StandardGapParams::cutoff3)
        .def_readwrite("cutoff3_width", &utils::StandardGapParams::cutoff3_width)
        .def_readwrite("n_sparse3", &utils::StandardGapParams::n_sparse3)
        .def("__repr__", [](const utils::StandardGapParams& p) {
            std::ostringstream oss;
            oss << "<StandardGapParams cutoff2=" << p.cutoff2
                << ", n_sparse2=" << p.n_sparse2
                << ", eam_n_sparse=" << p.eam_n_sparse
                << ", cutoff3=" << p.cutoff3
                << ", n_sparse3=" << p.n_sparse3
                << ", seed=" << p.seed << ">";
            return oss.str();
        });

    // =========================================================================
    // Standard Gap Fit
    // =========================================================================
    m.def("standard_gap_fit", [](const std::vector<Atoms>& training_data, const utils::StandardGapParams& params) {
        auto pot = utils::standardGapFit(training_data, params);
        return pot;
    }, py::arg("training_data"), py::arg("params") = utils::StandardGapParams(), py::call_guard<py::gil_scoped_release>(),
       "Fit a standard GAP potential from training data");

    // =========================================================================
    // Potential loader functions
    // =========================================================================
    m.def("load_potential", [](const std::string& path) -> std::shared_ptr<Potential> {
        auto val = loadPotential(path);
        return std::shared_ptr<Potential>(val.release());
    }, py::arg("path"), "Load a Potential from a single file path or stem");

    m.def("load_potential", [](const std::vector<std::string>& paths) -> std::shared_ptr<Potential> {
        auto val = loadPotential(paths);
        return std::shared_ptr<Potential>(val.release());
    }, py::arg("paths"), "Load a Potential from multiple file paths (e.g. .tabgap.h5 and .eam.fs)");

    m.def("read_atoms", [](const std::string& filename) {
        return Atoms::readAtoms(filename);
    }, py::arg("filename"), "Read XYZ dataset into a list of Atoms");

    m.def("write_atoms", [](const std::vector<Atoms>& frames, const std::string& filename) {
        if (frames.empty()) return;
        frames[0].write(filename); // first frame overwrites
        for (size_t i = 1; i < frames.size(); ++i) {
            std::ofstream out(filename, std::ios::app);
            frames[i].write(out);
        }
    }, py::arg("frames"), py::arg("filename"), "Write a list of Atoms frames to an XYZ file");

    // =========================================================================
    // Legacy / direct ASE calculator
    // =========================================================================
    py::class_<PyJGAPCalculator>(m, "PyJGAPCalculator")
        .def(py::init<const std::vector<std::string>&>(), py::arg("paths"))
        .def("calculate", &PyJGAPCalculator::calculate);
}
