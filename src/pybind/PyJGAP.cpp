#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
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
#include "jgap/core/atomic/species/composition/Species2Sorted.hpp"
#include "jgap/core/atomic/species/composition/Species3AtomicSorted.hpp"
#include "jgap/core/fit/gap/regularization/PerConfigTypeRegularizationRules.hpp"
#include "jgap/core/fit/gap/regularization/PerConfigTypeSigmas.hpp"
#include "jgap/core/fit/gap/regularization/Regularization.hpp"
#include "jgap/core/fit/gap/regularization/RegularizationRules.hpp"
#include "jgap/core/fit/gap/regularization/ScaledRegularizationRules.hpp"
#include "jgap/core/fit/gap/regularization/SimpleRegularizationRules.hpp"
#include "jgap/core/potentials/Cutoffs.hpp"
#include "jgap/core/potentials/Potential.hpp"
#include "jgap/impl/transform/nbody/3b/Distances3bTransformation.hpp"
#include "jgap/core/io/xyz/XYZData.hpp"
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
        const double* ptr = static_cast<double*>(buf.ptr);
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
        const double* ptr = static_cast<double*>(buf.ptr);
        if (buf.ndim == 1 && buf.shape[0] == 6) {
            return Virials{ptr[0], ptr[5], ptr[4], ptr[1], ptr[3], ptr[2]};
        } else if (buf.ndim == 2 && buf.shape[0] == 3 && buf.shape[1] == 3) {
            return Virials{ptr[0], ptr[1], ptr[2], ptr[4], ptr[5], ptr[8]};
        }
        throw std::invalid_argument("Expected (6,) Voigt array or (3, 3) matrix for Virials");
    }

    inline std::optional<Species> parseSpeciesOpt(const py::object& obj) {
        if (obj.is_none()) return std::nullopt;
        if (py::isinstance<Species>(obj)) return obj.cast<Species>();
        if (py::isinstance<py::str>(obj)) return Species(obj.cast<std::string>());
        throw std::invalid_argument("Expected Species, str, or None");
    }

    inline std::optional<Species2Sorted> parseSpecies2Opt(const py::object& obj) {
        if (obj.is_none()) return std::nullopt;
        if (py::isinstance<Species2Sorted>(obj)) return obj.cast<Species2Sorted>();
        if (py::isinstance<py::str>(obj)) {
            std::string s = obj.cast<std::string>();
            std::replace(s.begin(), s.end(), '-', ',');
            return Species2Sorted(s);
        }
        if (py::isinstance<py::tuple>(obj) || py::isinstance<py::list>(obj)) {
            auto seq = obj.cast<py::sequence>();
            if (seq.size() != 2) throw std::invalid_argument("Expected 2 species elements");
            Species s1 = py::isinstance<Species>(seq[0]) ? seq[0].cast<Species>() : Species(seq[0].cast<std::string>());
            Species s2 = py::isinstance<Species>(seq[1]) ? seq[1].cast<Species>() : Species(seq[1].cast<std::string>());
            return Species2Sorted(s1, s2);
        }
        throw std::invalid_argument("Expected Species2Sorted, str, 2-tuple, or None");
    }

    inline std::optional<Species3AtomicSorted> parseSpecies3Opt(const py::object& obj) {
        if (obj.is_none()) return std::nullopt;
        if (py::isinstance<Species3AtomicSorted>(obj)) return obj.cast<Species3AtomicSorted>();
        if (py::isinstance<py::str>(obj)) {
            std::string s = obj.cast<std::string>();
            if (s.find('|') == std::string::npos) {
                std::replace(s.begin(), s.end(), '-', ',');
                auto comma1 = s.find(',');
                if (comma1 != std::string::npos) {
                    std::string root = s.substr(0, comma1);
                    std::string rest = s.substr(comma1 + 1);
                    s = root + "|" + rest;
                }
            } else {
                std::replace(s.begin(), s.end(), '-', ',');
            }
            return Species3AtomicSorted(s);
        }
        if (py::isinstance<py::tuple>(obj) || py::isinstance<py::list>(obj)) {
            auto seq = obj.cast<py::sequence>();
            if (seq.size() != 3) throw std::invalid_argument("Expected 3 species elements (root, node1, node2)");
            Species root = py::isinstance<Species>(seq[0]) ? seq[0].cast<Species>() : Species(seq[0].cast<std::string>());
            Species s1 = py::isinstance<Species>(seq[1]) ? seq[1].cast<Species>() : Species(seq[1].cast<std::string>());
            Species s2 = py::isinstance<Species>(seq[2]) ? seq[2].cast<Species>() : Species(seq[2].cast<std::string>());
            return Species3AtomicSorted(root, s1, s2);
        }
        throw std::invalid_argument("Expected Species3AtomicSorted, str, 3-tuple, or None");
    }

    inline std::vector<utils::StandardGap2bParams> parse2bParamsList(const py::object& obj) {
        std::vector<utils::StandardGap2bParams> res;
        if (obj.is_none()) return res;
        for (const auto& item : obj) {
            res.push_back(item.cast<utils::StandardGap2bParams>());
        }
        return res;
    }

    inline std::vector<utils::StandardGapEamParams> parseEamParamsList(const py::object& obj) {
        std::vector<utils::StandardGapEamParams> res;
        if (obj.is_none()) return res;
        for (const auto& item : obj) {
            res.push_back(item.cast<utils::StandardGapEamParams>());
        }
        return res;
    }

    inline std::vector<utils::StandardGap3bParams> parse3bParamsList(const py::object& obj) {
        std::vector<utils::StandardGap3bParams> res;
        if (obj.is_none()) return res;
        for (const auto& item : obj) {
            res.push_back(item.cast<utils::StandardGap3bParams>());
        }
        return res;
    }

} // anonymous namespace

PYBIND11_MODULE(_jgap, m) {
    m.doc() = "jgap Python C++ bindings";

    // =========================================================================
    // Vector3
    // =========================================================================
    py::class_<Vector3>(m, "Vector3")
        .def(py::init<>())
        .def(py::init<double, double, double>(), py::arg("x"), py::arg("y"), py::arg("z"))
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
        .def("__mul__", [](const Vector3& a, double s) { return a * s; })
        .def("__rmul__", [](const Vector3& a, double s) { return a * s; })
        .def("__truediv__", [](const Vector3& a, double s) { return a / s; })
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
    // Species2Sorted & Species3AtomicSorted
    // =========================================================================
    py::class_<Species2Sorted>(m, "Species2Sorted")
        .def(py::init<Species, Species>(), py::arg("s1"), py::arg("s2"))
        .def(py::init<const std::string&>(), py::arg("encoded"))
        .def_property_readonly("nodes", [](const Species2Sorted& s) {
            return py::make_tuple(s.nodes[0], s.nodes[1]);
        })
        .def("to_string", &Species2Sorted::toString)
        .def("__str__", &Species2Sorted::toString)
        .def("__repr__", [](const Species2Sorted& s) {
            return "<Species2Sorted '" + s.toString() + "'>";
        })
        .def("__eq__", &Species2Sorted::operator==)
        .def("__lt__", &Species2Sorted::operator<)
        .def("__hash__", [](const Species2Sorted& s) {
            return std::hash<std::string>{}(s.toString());
        });

    py::class_<Species3AtomicSorted>(m, "Species3AtomicSorted")
        .def(py::init<Species, Species, Species>(), py::arg("root"), py::arg("s1"), py::arg("s2"))
        .def(py::init<const std::string&>(), py::arg("encoded"))
        .def_property_readonly("root", [](const Species3AtomicSorted& s) { return s.root; })
        .def_property_readonly("nodes", [](const Species3AtomicSorted& s) {
            return py::make_tuple(s.nodes[0], s.nodes[1]);
        })
        .def("to_string", &Species3AtomicSorted::toString)
        .def("__str__", &Species3AtomicSorted::toString)
        .def("__repr__", [](const Species3AtomicSorted& s) {
            return "<Species3AtomicSorted '" + s.toString() + "'>";
        })
        .def("__eq__", &Species3AtomicSorted::operator==)
        .def("__lt__", &Species3AtomicSorted::operator<)
        .def("__hash__", [](const Species3AtomicSorted& s) {
            return std::hash<std::string>{}(s.toString());
        });

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
        .def(py::init<double, double, double, double, double, double>(),
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
    // MainXYZPropertyNames
    // =========================================================================
    py::class_<MainXYZPropertyNames>(m, "MainXYZPropertyNames")
        .def(py::init<std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string>(),
             py::arg("positions") = "pos",
             py::arg("species") = "species",
             py::arg("forces") = "force",
             py::arg("virials") = "virial",
             py::arg("energy") = "energy",
             py::arg("lattice") = "Lattice",
             py::arg("pbc") = "pbc",
             py::arg("config_type") = "config_type")
        .def_readwrite("positions", &MainXYZPropertyNames::positions)
        .def_readwrite("species", &MainXYZPropertyNames::species)
        .def_readwrite("forces", &MainXYZPropertyNames::forces)
        .def_readwrite("virials", &MainXYZPropertyNames::virials)
        .def_readwrite("energy", &MainXYZPropertyNames::energy)
        .def_readwrite("lattice", &MainXYZPropertyNames::lattice)
        .def_readwrite("pbc", &MainXYZPropertyNames::pbc)
        .def_readwrite("config_type", &MainXYZPropertyNames::config_type)
        .def("__repr__", [](const MainXYZPropertyNames& p) {
            std::ostringstream oss;
            oss << "<MainXYZPropertyNames energy=\"" << p.energy
                << "\", forces=\"" << p.forces
                << "\", virials=\"" << p.virials << "\">";
            return oss.str();
        });

    auto readAtomsImpl = [](
        const std::string& filename,
        const std::optional<MainXYZPropertyNames>& prop_names,
        const std::optional<std::string>& positions,
        const std::optional<std::string>& species,
        const std::optional<std::string>& forces,
        const std::optional<std::string>& force,
        const std::optional<std::string>& virials,
        const std::optional<std::string>& virial,
        const std::optional<std::string>& energy,
        const std::optional<std::string>& lattice,
        const std::optional<std::string>& pbc,
        const std::optional<std::string>& config_type
    ) {
        MainXYZPropertyNames props = prop_names.value_or(MainXYZPropertyNames{});
        if (positions) props.positions = *positions;
        if (species) props.species = *species;
        if (forces) props.forces = *forces;
        if (force) props.forces = *force;
        if (virials) props.virials = *virials;
        if (virial) props.virials = *virial;
        if (energy) props.energy = *energy;
        if (lattice) props.lattice = *lattice;
        if (pbc) props.pbc = *pbc;
        if (config_type) props.config_type = *config_type;
        return Atoms::readAtoms(filename, props);
    };

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

        .def_property("energy", &Atoms::getEnergy, [](Atoms& a, std::optional<double> e) {
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

        .def_static("read_atoms", readAtomsImpl,
            py::arg("filename"),
            py::arg("prop_names") = std::nullopt,
            py::arg("positions") = std::nullopt,
            py::arg("species") = std::nullopt,
            py::arg("forces") = std::nullopt,
            py::arg("force") = std::nullopt,
            py::arg("virials") = std::nullopt,
            py::arg("virial") = std::nullopt,
            py::arg("energy") = std::nullopt,
            py::arg("lattice") = std::nullopt,
            py::arg("pbc") = std::nullopt,
            py::arg("config_type") = std::nullopt,
            py::call_guard<py::gil_scoped_release>(),
            "Read XYZ dataset into a list of Atoms, optionally specifying custom property names")

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
        .def(py::init<double>(), py::arg("energy"))
        .def(py::init<double, double, double>(), py::arg("energy"), py::arg("force"), py::arg("virials"))
        .def(py::init<double, double, double, double>(), py::arg("energy"), py::arg("force"), py::arg("virials_iso"), py::arg("virials_aniso"))
        .def(py::init<double, Vector3, Virials>(), py::arg("energy"), py::arg("force"), py::arg("virials"))
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
        .def(py::init<double, double, double, double, double, double>(),
             py::arg("energy_sigma_per_atom") = 0.001,
             py::arg("force_component_sigma") = 0.05,
             py::arg("virials_iso_sigma_per_atom") = 0.1,
             py::arg("virials_aniso_sigmas_per_atom") = 0.02,
             py::arg("liquid_multiplier") = 5.0,
             py::arg("short_range_multiplier") = 5.0)
        .def("__repr__", [](const SimpleRegularizationRules&) {
            return "<SimpleRegularizationRules>";
        });

    py::class_<ScaledRegularizationRules, RegularizationRules, std::shared_ptr<ScaledRegularizationRules>>(m, "ScaledRegularizationRules")
        .def(py::init<std::shared_ptr<RegularizationRules>, double, double>(),
             py::arg("base_rules"),
             py::arg("force_scale") = 1.0,
             py::arg("min_scale") = 1.0)
        .def(py::init<PerConfigTypeSigmas, double, double>(),
             py::arg("base_sigmas"),
             py::arg("force_scale") = 1.0,
             py::arg("min_scale") = 1.0)
        .def(py::init<double, double, double, double, double, double>(),
             py::arg("energy_sigma_per_atom") = 0.001,
             py::arg("force_component_sigma") = 0.05,
             py::arg("virials_iso_sigma_per_atom") = 0.1,
             py::arg("virials_aniso_sigmas_per_atom") = 0.02,
             py::arg("force_scale") = 1.0,
             py::arg("min_scale") = 1.0)
        .def_property_readonly("force_scale", &ScaledRegularizationRules::getForceScale)
        .def_property_readonly("min_scale", &ScaledRegularizationRules::getMinScale)
        .def("__repr__", [](const ScaledRegularizationRules& r) {
            return "<ScaledRegularizationRules force_scale=" + std::to_string(r.getForceScale()) + ">";
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
    // 3B Transformation Enum
    // =========================================================================
    py::enum_<utils::ThreeBodyTransformationType>(m, "ThreeBodyTransformationType")
        .value("Angle", utils::ThreeBodyTransformationType::Angle)
        .value("Distances", utils::ThreeBodyTransformationType::Distances)
        .export_values();

    // =========================================================================
    // StandardGap 2b, Eam, 3b Params
    // =========================================================================
    py::class_<utils::StandardGap2bParams>(m, "StandardGap2bParams")
        .def(py::init([](py::object species, double cutoff, double cutoff_width, size_t n_sparse, double energy_scale, double length_scale) {
            utils::StandardGap2bParams p;
            p.species = parseSpecies2Opt(species);
            p.cutoff = cutoff;
            p.cutoff_width = cutoff_width;
            p.n_sparse = n_sparse;
            p.energy_scale = energy_scale;
            p.length_scale = length_scale;
            return p;
        }),
        py::arg("species") = py::none(),
        py::arg("cutoff") = 4.5,
        py::arg("cutoff_width") = 1.0,
        py::arg("n_sparse") = 20,
        py::arg("energy_scale") = 10.0,
        py::arg("length_scale") = 1.0)
        .def_property("species",
            [](const utils::StandardGap2bParams& p) -> py::object {
                if (p.species) return py::cast(*p.species);
                return py::none();
            },
            [](utils::StandardGap2bParams& p, py::object val) {
                p.species = parseSpecies2Opt(val);
            })
        .def_readwrite("cutoff", &utils::StandardGap2bParams::cutoff)
        .def_readwrite("cutoff_width", &utils::StandardGap2bParams::cutoff_width)
        .def_readwrite("n_sparse", &utils::StandardGap2bParams::n_sparse)
        .def_readwrite("energy_scale", &utils::StandardGap2bParams::energy_scale)
        .def_readwrite("length_scale", &utils::StandardGap2bParams::length_scale)
        .def("__repr__", [](const utils::StandardGap2bParams& p) {
            std::ostringstream oss;
            oss << "<StandardGap2bParams species=" << (p.species ? p.species->toString() : "None")
                << ", cutoff=" << p.cutoff
                << ", cutoff_width=" << p.cutoff_width
                << ", n_sparse=" << p.n_sparse
                << ", energy_scale=" << p.energy_scale
                << ", length_scale=" << p.length_scale << ">";
            return oss.str();
        })
        .def("__eq__", &utils::StandardGap2bParams::operator==);

    py::class_<utils::StandardGapEamParams>(m, "StandardGapEamParams")
        .def(py::init([](py::object species, EamMode eam_mode, utils::EamPairFunctionType eam_pair_function, double cutoff, size_t n_sparse, double min_density, double energy_scale, double length_scale) {
            utils::StandardGapEamParams p;
            p.species = parseSpeciesOpt(species);
            p.eam_mode = eam_mode;
            p.eam_pair_function = eam_pair_function;
            p.cutoff = cutoff;
            p.n_sparse = n_sparse;
            p.min_density = min_density;
            p.energy_scale = energy_scale;
            p.length_scale = length_scale;
            return p;
        }),
        py::arg("species") = py::none(),
        py::arg("eam_mode") = EamMode::Blind,
        py::arg("eam_pair_function") = utils::EamPairFunctionType::FSGen3,
        py::arg("cutoff") = 4.5,
        py::arg("n_sparse") = 20,
        py::arg("min_density") = 0.05,
        py::arg("energy_scale") = 1.0,
        py::arg("length_scale") = 1.0)
        .def_property("species",
            [](const utils::StandardGapEamParams& p) -> py::object {
                if (p.species) return py::cast(*p.species);
                return py::none();
            },
            [](utils::StandardGapEamParams& p, py::object val) {
                p.species = parseSpeciesOpt(val);
            })
        .def_readwrite("eam_mode", &utils::StandardGapEamParams::eam_mode)
        .def_readwrite("eam_pair_function", &utils::StandardGapEamParams::eam_pair_function)
        .def_readwrite("cutoff", &utils::StandardGapEamParams::cutoff)
        .def_readwrite("n_sparse", &utils::StandardGapEamParams::n_sparse)
        .def_readwrite("min_density", &utils::StandardGapEamParams::min_density)
        .def_readwrite("energy_scale", &utils::StandardGapEamParams::energy_scale)
        .def_readwrite("length_scale", &utils::StandardGapEamParams::length_scale)
        .def_property("density_scale",
            [](const utils::StandardGapEamParams& p) { return p.length_scale; },
            [](utils::StandardGapEamParams& p, double v) { p.length_scale = v; })
        .def("__repr__", [](const utils::StandardGapEamParams& p) {
            std::ostringstream oss;
            oss << "<StandardGapEamParams species=" << (p.species ? p.species->symbol() : "None")
                << ", cutoff=" << p.cutoff
                << ", n_sparse=" << p.n_sparse
                << ", min_density=" << p.min_density
                << ", energy_scale=" << p.energy_scale
                << ", length_scale=" << p.length_scale << ">";
            return oss.str();
        })
        .def("__eq__", &utils::StandardGapEamParams::operator==);

    py::class_<utils::StandardGap3bParams>(m, "StandardGap3bParams")
        .def(py::init([](py::object species, utils::ThreeBodyTransformationType transformation_type, double cutoff, double cutoff_width, size_t n_sparse, double energy_scale, py::object length_scales, py::object length_scale) {
            utils::StandardGap3bParams p;
            p.species = parseSpecies3Opt(species);
            p.transformation_type = transformation_type;
            p.cutoff = cutoff;
            p.cutoff_width = cutoff_width;
            p.n_sparse = n_sparse;
            p.energy_scale = energy_scale;
            if (!length_scale.is_none()) {
                double val = length_scale.cast<double>();
                p.length_scales = {val, val, val};
            }
            if (!length_scales.is_none()) {
                if (py::isinstance<py::float_>(length_scales) || py::isinstance<py::int_>(length_scales)) {
                    double val = length_scales.cast<double>();
                    p.length_scales = {val, val, val};
                } else {
                    auto vec = length_scales.cast<std::vector<double>>();
                    if (vec.size() != 3) throw std::invalid_argument("Expected length_scales of size 3");
                    p.length_scales = {vec[0], vec[1], vec[2]};
                }
            }
            return p;
        }),
        py::arg("species") = py::none(),
        py::arg("transformation_type") = utils::ThreeBodyTransformationType::Angle,
        py::arg("cutoff") = 3.7,
        py::arg("cutoff_width") = 0.6,
        py::arg("n_sparse") = 500,
        py::arg("energy_scale") = 1.0,
        py::arg("length_scales") = py::none(),
        py::arg("length_scale") = py::none())
        .def_property("species",
            [](const utils::StandardGap3bParams& p) -> py::object {
                if (p.species) return py::cast(*p.species);
                return py::none();
            },
            [](utils::StandardGap3bParams& p, py::object val) {
                p.species = parseSpecies3Opt(val);
            })
        .def_readwrite("transformation_type", &utils::StandardGap3bParams::transformation_type)
        .def_readwrite("cutoff", &utils::StandardGap3bParams::cutoff)
        .def_readwrite("cutoff_width", &utils::StandardGap3bParams::cutoff_width)
        .def_readwrite("n_sparse", &utils::StandardGap3bParams::n_sparse)
        .def_readwrite("energy_scale", &utils::StandardGap3bParams::energy_scale)
        .def_property("length_scales",
            [](const utils::StandardGap3bParams& p) { return p.length_scales; },
            [](utils::StandardGap3bParams& p, const std::array<double, 3>& v) { p.length_scales = v; })
        .def_property("length_scale",
            [](const utils::StandardGap3bParams& p) { return p.length_scales[0]; },
            [](utils::StandardGap3bParams& p, double v) { p.length_scales = {v, v, v}; })
        .def("__repr__", [](const utils::StandardGap3bParams& p) {
            std::ostringstream oss;
            oss << "<StandardGap3bParams species=" << (p.species ? p.species->toString() : "None")
                << ", transformation_type=" << (p.transformation_type == utils::ThreeBodyTransformationType::Angle ? "Angle" : "Distances")
                << ", cutoff=" << p.cutoff
                << ", cutoff_width=" << p.cutoff_width
                << ", n_sparse=" << p.n_sparse
                << ", energy_scale=" << p.energy_scale
                << ", length_scales=[" << p.length_scales[0] << ", " << p.length_scales[1] << ", " << p.length_scales[2] << "]>";
            return oss.str();
        })
        .def("__eq__", &utils::StandardGap3bParams::operator==);

    // =========================================================================
    // StandardGapParams
    // =========================================================================
    py::class_<utils::StandardGapParams>(m, "StandardGapParams")
        .def(py::init([](
            size_t seed,
            std::optional<std::string> screened_coulomb_dataset_file,
            double approx_ram_limit_gb,
            std::optional<utils::StandardGap2bParams> default_2b,
            std::optional<utils::StandardGapEamParams> default_eam,
            std::optional<utils::StandardGap3bParams> default_3b,
            py::object species_2b,
            py::object species_eam,
            py::object species_3b,
            // Legacy kwargs for backward compatibility
            std::optional<double> cutoff2,
            std::optional<double> cutoff2_width,
            std::optional<size_t> n_sparse2,
            std::optional<EamMode> eam_mode,
            std::optional<utils::EamPairFunctionType> eam_pair_function,
            std::optional<size_t> eam_n_sparse,
            std::optional<double> eam_min_density,
            std::optional<double> cutoff3,
            std::optional<double> cutoff3_width,
            std::optional<size_t> n_sparse3,
            std::optional<double> energy_scale_2b,
            std::optional<double> length_scale_2b,
            std::optional<double> energy_scale_eam,
            std::optional<double> length_scale_eam,
            std::optional<double> energy_scale_3b,
            py::object length_scales_3b,
            std::optional<double> length_scale_3b
        ) {
            utils::StandardGapParams p;
            p.seed = seed;
            p.screened_coulomb_dataset_file = screened_coulomb_dataset_file;
            p.approx_ram_limit_gb = approx_ram_limit_gb;
            p.default_2b = default_2b;
            p.default_eam = default_eam;
            p.default_3b = default_3b;
            p.species_2b = parse2bParamsList(species_2b);
            p.species_eam = parseEamParamsList(species_eam);
            p.species_3b = parse3bParamsList(species_3b);

            // Apply legacy kwargs if default was retained
            if (p.default_2b) {
                if (cutoff2.has_value()) p.default_2b->cutoff = *cutoff2;
                if (cutoff2_width.has_value()) p.default_2b->cutoff_width = *cutoff2_width;
                if (n_sparse2.has_value()) p.default_2b->n_sparse = *n_sparse2;
            }
            if (p.default_eam) {
                if (eam_mode.has_value()) p.default_eam->eam_mode = *eam_mode;
                if (eam_pair_function.has_value()) p.default_eam->eam_pair_function = *eam_pair_function;
                if (eam_n_sparse.has_value()) p.default_eam->n_sparse = *eam_n_sparse;
                if (eam_min_density.has_value()) p.default_eam->min_density = *eam_min_density;
            }
            if (p.default_2b) {
                if (energy_scale_2b.has_value()) p.default_2b->energy_scale = *energy_scale_2b;
                if (length_scale_2b.has_value()) p.default_2b->length_scale = *length_scale_2b;
            }
            if (p.default_eam) {
                if (energy_scale_eam.has_value()) p.default_eam->energy_scale = *energy_scale_eam;
                if (length_scale_eam.has_value()) p.default_eam->length_scale = *length_scale_eam;
            }
            if (p.default_3b) {
                if (cutoff3.has_value()) p.default_3b->cutoff = *cutoff3;
                if (cutoff3_width.has_value()) p.default_3b->cutoff_width = *cutoff3_width;
                if (n_sparse3.has_value()) p.default_3b->n_sparse = *n_sparse3;
                if (energy_scale_3b.has_value()) p.default_3b->energy_scale = *energy_scale_3b;
                if (length_scale_3b.has_value()) p.default_3b->length_scales = {*length_scale_3b, *length_scale_3b, *length_scale_3b};
                if (!length_scales_3b.is_none()) {
                    if (py::isinstance<py::float_>(length_scales_3b) || py::isinstance<py::int_>(length_scales_3b)) {
                        double v = length_scales_3b.cast<double>();
                        p.default_3b->length_scales = {v, v, v};
                    } else {
                        auto vec = length_scales_3b.cast<std::vector<double>>();
                        if (vec.size() == 3) p.default_3b->length_scales = {vec[0], vec[1], vec[2]};
                    }
                }
            }
            return p;
        }),
        py::arg("seed") = 42,
        py::arg("screened_coulomb_dataset_file") = std::nullopt,
        py::arg("approx_ram_limit_gb") = 4.0,
        py::arg("default_2b") = utils::StandardGap2bParams{},
        py::arg("default_eam") = utils::StandardGapEamParams{},
        py::arg("default_3b") = utils::StandardGap3bParams{},
        py::arg("species_2b") = py::none(),
        py::arg("species_eam") = py::none(),
        py::arg("species_3b") = py::none(),
        py::arg("cutoff2") = std::nullopt,
        py::arg("cutoff2_width") = std::nullopt,
        py::arg("n_sparse2") = std::nullopt,
        py::arg("eam_mode") = std::nullopt,
        py::arg("eam_pair_function") = std::nullopt,
        py::arg("eam_n_sparse") = std::nullopt,
        py::arg("eam_min_density") = std::nullopt,
        py::arg("cutoff3") = std::nullopt,
        py::arg("cutoff3_width") = std::nullopt,
        py::arg("n_sparse3") = std::nullopt,
        py::arg("energy_scale_2b") = std::nullopt,
        py::arg("length_scale_2b") = std::nullopt,
        py::arg("energy_scale_eam") = std::nullopt,
        py::arg("length_scale_eam") = std::nullopt,
        py::arg("energy_scale_3b") = std::nullopt,
        py::arg("length_scales_3b") = py::none(),
        py::arg("length_scale_3b") = std::nullopt)
        .def_readwrite("seed", &utils::StandardGapParams::seed)
        .def_readwrite("screened_coulomb_dataset_file", &utils::StandardGapParams::screened_coulomb_dataset_file)
        .def_readwrite("approx_ram_limit_gb", &utils::StandardGapParams::approx_ram_limit_gb)
        .def_readwrite("default_2b", &utils::StandardGapParams::default_2b)
        .def_readwrite("default_eam", &utils::StandardGapParams::default_eam)
        .def_readwrite("default_3b", &utils::StandardGapParams::default_3b)
        .def_property(
            "species_2b",
            [](const utils::StandardGapParams& p) { return p.species_2b; },
            [](utils::StandardGapParams& p, const py::object& obj) { p.species_2b = parse2bParamsList(obj); })
        .def_property(
            "species_eam",
            [](const utils::StandardGapParams& p) { return p.species_eam; },
            [](utils::StandardGapParams& p, const py::object& obj) { p.species_eam = parseEamParamsList(obj); })
        .def_property(
            "species_3b",
            [](const utils::StandardGapParams& p) { return p.species_3b; },
            [](utils::StandardGapParams& p, const py::object& obj) { p.species_3b = parse3bParamsList(obj); })
        .def("add_species_2b", [](utils::StandardGapParams& p, const utils::StandardGap2bParams& item) {
            p.species_2b.push_back(item);
        })
        .def("add_species_eam", [](utils::StandardGapParams& p, const utils::StandardGapEamParams& item) {
            p.species_eam.push_back(item);
        })
        .def("add_species_3b", [](utils::StandardGapParams& p, const utils::StandardGap3bParams& item) {
            p.species_3b.push_back(item);
        })
        // Backwards compatibility properties forwarding to default_2b / default_eam / default_3b
        .def_property("cutoff2",
            [](const utils::StandardGapParams& p) { return p.default_2b ? p.default_2b->cutoff : 0.0; },
            [](utils::StandardGapParams& p, double v) { if (p.default_2b) p.default_2b->cutoff = v; })
        .def_property("cutoff2_width",
            [](const utils::StandardGapParams& p) { return p.default_2b ? p.default_2b->cutoff_width : 0.0; },
            [](utils::StandardGapParams& p, double v) { if (p.default_2b) p.default_2b->cutoff_width = v; })
        .def_property("n_sparse2",
            [](const utils::StandardGapParams& p) { return p.default_2b ? p.default_2b->n_sparse : 0; },
            [](utils::StandardGapParams& p, size_t v) { if (p.default_2b) p.default_2b->n_sparse = v; })
        .def_property("eam_mode",
            [](const utils::StandardGapParams& p) { return p.default_eam ? p.default_eam->eam_mode : EamMode::Blind; },
            [](utils::StandardGapParams& p, EamMode v) { if (p.default_eam) p.default_eam->eam_mode = v; })
        .def_property("eam_pair_function",
            [](const utils::StandardGapParams& p) { return p.default_eam ? p.default_eam->eam_pair_function : utils::EamPairFunctionType::FSGen3; },
            [](utils::StandardGapParams& p, utils::EamPairFunctionType v) { if (p.default_eam) p.default_eam->eam_pair_function = v; })
        .def_property("eam_n_sparse",
            [](const utils::StandardGapParams& p) { return p.default_eam ? p.default_eam->n_sparse : 0; },
            [](utils::StandardGapParams& p, size_t v) { if (p.default_eam) p.default_eam->n_sparse = v; })
        .def_property("eam_min_density",
            [](const utils::StandardGapParams& p) { return p.default_eam ? p.default_eam->min_density : 0.0; },
            [](utils::StandardGapParams& p, double v) { if (p.default_eam) p.default_eam->min_density = v; })
        .def_property("cutoff3",
            [](const utils::StandardGapParams& p) { return p.default_3b ? p.default_3b->cutoff : 0.0; },
            [](utils::StandardGapParams& p, double v) { if (p.default_3b) p.default_3b->cutoff = v; })
        .def_property("cutoff3_width",
            [](const utils::StandardGapParams& p) { return p.default_3b ? p.default_3b->cutoff_width : 0.0; },
            [](utils::StandardGapParams& p, double v) { if (p.default_3b) p.default_3b->cutoff_width = v; })
        .def_property("n_sparse3",
            [](const utils::StandardGapParams& p) { return p.default_3b ? p.default_3b->n_sparse : 0; },
            [](utils::StandardGapParams& p, size_t v) { if (p.default_3b) p.default_3b->n_sparse = v; })
        .def_property("energy_scale_2b",
            [](const utils::StandardGapParams& p) { return p.default_2b ? p.default_2b->energy_scale : 10.0; },
            [](utils::StandardGapParams& p, double v) { if (p.default_2b) p.default_2b->energy_scale = v; })
        .def_property("length_scale_2b",
            [](const utils::StandardGapParams& p) { return p.default_2b ? p.default_2b->length_scale : 1.0; },
            [](utils::StandardGapParams& p, double v) { if (p.default_2b) p.default_2b->length_scale = v; })
        .def_property("energy_scale_eam",
            [](const utils::StandardGapParams& p) { return p.default_eam ? p.default_eam->energy_scale : 1.0; },
            [](utils::StandardGapParams& p, double v) { if (p.default_eam) p.default_eam->energy_scale = v; })
        .def_property("length_scale_eam",
            [](const utils::StandardGapParams& p) { return p.default_eam ? p.default_eam->length_scale : 1.0; },
            [](utils::StandardGapParams& p, double v) { if (p.default_eam) p.default_eam->length_scale = v; })
        .def_property("energy_scale_3b",
            [](const utils::StandardGapParams& p) { return p.default_3b ? p.default_3b->energy_scale : 1.0; },
            [](utils::StandardGapParams& p, double v) { if (p.default_3b) p.default_3b->energy_scale = v; })
        .def_property("length_scale_3b",
            [](const utils::StandardGapParams& p) { return p.default_3b ? p.default_3b->length_scales[0] : 1.0; },
            [](utils::StandardGapParams& p, double v) { if (p.default_3b) p.default_3b->length_scales = {v, v, v}; })
        .def_property("length_scales_3b",
            [](const utils::StandardGapParams& p) -> py::object {
                if (p.default_3b) return py::cast(p.default_3b->length_scales);
                return py::none();
            },
            [](utils::StandardGapParams& p, const std::array<double, 3>& v) { if (p.default_3b) p.default_3b->length_scales = v; })
        .def("__repr__", [](const utils::StandardGapParams& p) {
            std::ostringstream oss;
            oss << "<StandardGapParams has_default_2b=" << p.default_2b.has_value()
                << ", has_default_eam=" << p.default_eam.has_value()
                << ", has_default_3b=" << p.default_3b.has_value()
                << ", species_2b=" << p.species_2b.size()
                << ", species_eam=" << p.species_eam.size()
                << ", species_3b=" << p.species_3b.size()
                << ", approx_ram_limit_gb=" << p.approx_ram_limit_gb
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
        .def(py::init<double, double, size_t, std::array<size_t, 3>>(),
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

    m.def("read_atoms", readAtomsImpl,
        py::arg("filename"),
        py::arg("prop_names") = std::nullopt,
        py::arg("positions") = std::nullopt,
        py::arg("species") = std::nullopt,
        py::arg("forces") = std::nullopt,
        py::arg("force") = std::nullopt,
        py::arg("virials") = std::nullopt,
        py::arg("virial") = std::nullopt,
        py::arg("energy") = std::nullopt,
        py::arg("lattice") = std::nullopt,
        py::arg("pbc") = std::nullopt,
        py::arg("config_type") = std::nullopt,
        py::call_guard<py::gil_scoped_release>(),
        "Read XYZ dataset into a list of Atoms, optionally specifying custom property names");

    m.def("write_atoms", [](const std::vector<Atoms>& frames, const std::string& filename) {
        if (frames.empty()) return;
        frames[0].write(filename); // first frame overwrites
        for (size_t i = 1; i < frames.size(); ++i) {
            std::ofstream out(filename, std::ios::app);
            frames[i].write(out);
        }
    }, py::arg("frames"), py::arg("filename"), "Write a list of Atoms frames to an XYZ file");
}
