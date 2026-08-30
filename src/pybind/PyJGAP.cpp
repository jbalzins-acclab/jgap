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
#include "jgap/core/fit/gap/regularization/PerConfigTypeRegularizationRules.hpp"
#include "jgap/core/fit/gap/regularization/PerConfigTypeSigmas.hpp"
#include "jgap/core/fit/gap/regularization/Regularization.hpp"
#include "jgap/core/fit/gap/regularization/RegularizationRules.hpp"
#include "jgap/core/fit/gap/regularization/SimpleRegularizationRules.hpp"
#include "jgap/core/potentials/Cutoffs.hpp"
#include "jgap/core/potentials/Potential.hpp"
#include "jgap/io/PotentialLoader.hpp"
#include "jgap/utils/gap/StandardGapFit.hpp"
#include "jgap/utils/gap/StandardGapParams.hpp"
#include "jgap/utils/gap/StandardTabulation.hpp"

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
    // Potential
    // =========================================================================
    py::class_<Potential, std::shared_ptr<Potential>> pot_class(m, "Potential");
    pot_class
        .def("calculate_energy", &Potential::calculateEnergy,
             py::arg("atoms"),
             py::call_guard<py::gil_scoped_release>(),
             "Calculate energy, forces, and virials for the given structure")
        .def("get_cutoffs", &Potential::getCutoffs)
        .def_static("load", [](const std::string& path) -> std::shared_ptr<Potential> {
            auto val = loadPotential(path);
            return std::shared_ptr<Potential>(val.release());
        }, py::arg("path"), "Load potential from file path or base stem");

    // =========================================================================
    // Regularization & Sigmas
    // =========================================================================
    py::class_<PerConfigTypeSigmas>(m, "PerConfigTypeSigmas")
        .def(py::init<Real>(), py::arg("energy"))
        .def(py::init<Real, Real, Real>(), py::arg("energy"), py::arg("force"), py::arg("virials"))
        .def(py::init<Real, Real, Real, Real>(), py::arg("energy"), py::arg("force"), py::arg("virials_iso"), py::arg("virials_aniso"))
        .def(py::init<Real, Vector3, Virials>(), py::arg("energy"), py::arg("force"), py::arg("virials"))
        .def_readwrite("energy", &PerConfigTypeSigmas::energy)
        .def_readwrite("force", &PerConfigTypeSigmas::force)
        .def_readwrite("virials", &PerConfigTypeSigmas::virials)
        .def("__repr__", [](const PerConfigTypeSigmas& s) {
            std::ostringstream oss;
            oss << "<PerConfigTypeSigmas energy=" << s.energy << ", force=" << s.force.x << ">";
            return oss.str();
        });

    py::class_<Regularization>(m, "Regularization")
        .def(py::init<>())
        .def_readwrite("energy", &Regularization::energy)
        .def_readwrite("virials", &Regularization::virials)
        .def_readwrite("forces", &Regularization::forces)
        .def("__repr__", [](const Regularization& r) {
            std::ostringstream oss;
            oss << "<Regularization has_energy=" << r.energy.has_value()
                << ", has_forces=" << r.forces.has_value()
                << ", has_virials=" << r.virials.has_value() << ">";
            return oss.str();
        });

    // =========================================================================
    // Regularization Rules
    // =========================================================================
    py::class_<RegularizationRules, std::shared_ptr<RegularizationRules>>(m, "RegularizationRules")
        .def("determine", &RegularizationRules::determine, py::arg("atoms"))
        .def("determine_for_all", &RegularizationRules::determineForAll, py::arg("structures"));

    py::class_<PerConfigTypeRegularizationRules, RegularizationRules, std::shared_ptr<PerConfigTypeRegularizationRules>>(m, "PerConfigTypeRegularizationRules")
        .def(py::init<PerConfigTypeSigmas, std::map<std::string, PerConfigTypeSigmas>, std::map<std::string, PerConfigTypeSigmas>>(),
             py::arg("default_sigmas"),
             py::arg("exact_config_type_sigmas") = std::map<std::string, PerConfigTypeSigmas>{},
             py::arg("config_type_contains_sigmas") = std::map<std::string, PerConfigTypeSigmas>{})
        .def(py::init<PerConfigTypeSigmas, const std::string&>(),
             py::arg("default_sigmas"),
             py::arg("config_string"))
        .def_property_readonly("defaults", &PerConfigTypeRegularizationRules::getDefaults)
        .def_property_readonly("exact_config_type_sigmas", &PerConfigTypeRegularizationRules::getExactConfigTypeSigmas)
        .def_property_readonly("config_type_contains_sigmas", &PerConfigTypeRegularizationRules::getConfigTypeContainsSigmas)
        .def("__repr__", [](const PerConfigTypeRegularizationRules&) {
            return "<PerConfigTypeRegularizationRules>";
        });

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
    // EAM Enums
    // =========================================================================
    py::enum_<utils::EamPairFunctionType>(m, "EamPairFunctionType")
        .value("FSGen2", utils::EamPairFunctionType::FSGen2)
        .value("FSGen3", utils::EamPairFunctionType::FSGen3)
        .value("Coscutoff", utils::EamPairFunctionType::Coscutoff)
        .value("Polycutoff", utils::EamPairFunctionType::Polycutoff)
        .export_values();

    py::enum_<EamMode>(m, "EamMode")
        .value("FSsym", EamMode::FSsym)
        .value("FSgen", EamMode::FSgen)
        .value("EAM", EamMode::EAM)
        .value("Blind", EamMode::Blind)
        .export_values();

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
            utils::EamPairFunctionType eam_pair_function,
            size_t eam_n_sparse,
            Real eam_min_density,
            Real cutoff3,
            Real cutoff3_width,
            size_t n_sparse3,
            std::optional<Real> approx_ram_limit_gb
        ) {
            utils::StandardGapParams p;
            p.seed = seed;
            p.screened_coulomb_dataset_file = screened_coulomb_dataset_file;
            p.cutoff2 = cutoff2;
            p.cutoff2_width = cutoff2_width;
            p.n_sparse2 = n_sparse2;
            p.eam_mode = eam_mode;
            p.eam_pair_function = eam_pair_function;
            p.eam_n_sparse = eam_n_sparse;
            p.eam_min_density = eam_min_density;
            p.cutoff3 = cutoff3;
            p.cutoff3_width = cutoff3_width;
            p.n_sparse3 = n_sparse3;
            p.approx_ram_limit_gb = approx_ram_limit_gb;
            return p;
        }),
        py::arg("seed") = 42,
        py::arg("screened_coulomb_dataset_file") = std::nullopt,
        py::arg("cutoff2") = 4.5,
        py::arg("cutoff2_width") = 1.0,
        py::arg("n_sparse2") = 20,
        py::arg("eam_mode") = EamMode::Blind,
        py::arg("eam_pair_function") = utils::EamPairFunctionType::FSGen3,
        py::arg("eam_n_sparse") = 20,
        py::arg("eam_min_density") = 0.05,
        py::arg("cutoff3") = 3.7,
        py::arg("cutoff3_width") = 0.6,
        py::arg("n_sparse3") = 500,
        py::arg("approx_ram_limit_gb") = std::nullopt)
        .def_readwrite("seed", &utils::StandardGapParams::seed)
        .def_readwrite("screened_coulomb_dataset_file", &utils::StandardGapParams::screened_coulomb_dataset_file)
        .def_readwrite("cutoff2", &utils::StandardGapParams::cutoff2)
        .def_readwrite("cutoff2_width", &utils::StandardGapParams::cutoff2_width)
        .def_readwrite("n_sparse2", &utils::StandardGapParams::n_sparse2)
        .def_readwrite("eam_mode", &utils::StandardGapParams::eam_mode)
        .def_readwrite("eam_pair_function", &utils::StandardGapParams::eam_pair_function)
        .def_readwrite("eam_n_sparse", &utils::StandardGapParams::eam_n_sparse)
        .def_readwrite("eam_min_density", &utils::StandardGapParams::eam_min_density)
        .def_readwrite("cutoff3", &utils::StandardGapParams::cutoff3)
        .def_readwrite("cutoff3_width", &utils::StandardGapParams::cutoff3_width)
        .def_readwrite("n_sparse3", &utils::StandardGapParams::n_sparse3)
        .def_readwrite("approx_ram_limit_gb", &utils::StandardGapParams::approx_ram_limit_gb)
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
    m.def("standard_gap_fit", [](
        const std::string& filename,
        const std::vector<Atoms>& training_data,
        const std::vector<Regularization>& sigmas,
        const utils::StandardGapParams& params
    ) {
        utils::standardGapFit(filename, training_data, sigmas, params);
    },
    py::arg("filename"),
    py::arg("training_data"),
    py::arg("sigmas"),
    py::arg("params") = utils::StandardGapParams(),
    py::call_guard<py::gil_scoped_release>(),
    "Fit a standard GAP potential from training data and write to file");

    // =========================================================================
    // Tabulation
    // =========================================================================
    py::class_<utils::StandardTabulationParams>(m, "StandardTabulationParams")
        .def(py::init<Real, Real, size_t, std::array<size_t, 3>>(),
             py::arg("r_min_3b") = 0.5,
             py::arg("max_eam_density") = 10.0,
             py::arg("n_grid_2b") = 5000,
             py::arg("n_grid_3b") = std::array<size_t, 3>{80, 80, 80})
        .def_readwrite("r_min_3b", &utils::StandardTabulationParams::r_min_3b)
        .def_readwrite("max_eam_density", &utils::StandardTabulationParams::max_eam_density)
        .def_readwrite("n_grid_2b", &utils::StandardTabulationParams::n_grid_2b)
        .def_readwrite("n_grid_3b", &utils::StandardTabulationParams::n_grid_3b)
        .def("__repr__", [](const utils::StandardTabulationParams& p) {
            std::ostringstream oss;
            oss << "<StandardTabulationParams r_min_3b=" << p.r_min_3b
                << ", max_eam_density=" << p.max_eam_density
                << ", n_grid_2b=" << p.n_grid_2b
                << ", n_grid_3b=[" << p.n_grid_3b[0] << "," << p.n_grid_3b[1] << "," << p.n_grid_3b[2] << "]>";
            return oss.str();
        });

    m.def("standard_tabulation", [](
        const std::string& pot_filename,
        const std::string& output_prefix,
        const utils::StandardTabulationParams& params
    ) {
        utils::standardTabulation(pot_filename, output_prefix, params);
    },
    py::arg("pot_filename"),
    py::arg("output_prefix"),
    py::arg("params") = utils::StandardTabulationParams(),
    py::call_guard<py::gil_scoped_release>(),
    "Tabulate a potential file and save .tabgap.h5 and .eam.fs files");

    m.def("standard_tabulation", [](
        const Potential& potential,
        const std::string& output_prefix,
        const utils::StandardTabulationParams& params
    ) {
        utils::standardTabulation(potential, output_prefix, params);
    },
    py::arg("potential"),
    py::arg("output_prefix"),
    py::arg("params") = utils::StandardTabulationParams(),
    py::call_guard<py::gil_scoped_release>(),
    "Tabulate a Potential object and save .tabgap.h5 and .eam.fs files");

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
}
