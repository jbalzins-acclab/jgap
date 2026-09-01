#ifndef JGAP_XYZDATA_HPP
#define JGAP_XYZDATA_HPP
#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <variant>
#include <vector>

#include "../../Vector3.hpp"
#include "jgap/core/Real.hpp"
#include "jgap/core/atomic/energy/Virials.hpp"
#include "jgap/core/atomic/geometry/Lattice.hpp"
#include "jgap/core/atomic/species/Species.hpp"

namespace jgap {
    using XYZInfoType = std::variant<std::string, int, Real, Vector3, Virials, Lattice, std::array<bool, 3> >;
    using XYZArrayType = std::variant<
        std::vector<int>,
        std::vector<Real>,
        std::vector<Vector3>,
        std::vector<std::string>,
        std::vector<Species> >;

    /// Labels for the common ext-xyz property/array names that @ref XYZData::read and @ref Atoms use.
    ///
    /// Defaults follow standard ASE naming for convenience,
    /// however, virials require as there seems to be no standard label.
    /// Otherwise, one might need to check "force" <-> "forces" ambiguity and "Lattice" casing.
    struct MainXYZPropertyNames {
        std::string positions = "pos";
        std::string species = "species";
        std::string forces = "force";
        std::string virials = "virial";
        std::string energy = "energy";
        std::string lattice = "Lattice";
        std::string pbc = "pbc";
        std::string config_type = "config_type";
    };

    /// Extended .xyz file data.
    /// Doesn't assume what exactly is known pre-atom,
    /// rather stores all the header properties,
    /// in some of the @ref XYZInfoType types,
    /// and a bunch of named arrays as @ref XYZArrayType.
    ///
    /// @see @ref XYZData::read for specifics.
    ///
    /// @note inspired by ASE.
    class XYZData {
    public:
        static std::vector<XYZData> read(const std::string& filename, const MainXYZPropertyNames& main_props = {});

        /// Tries to read the file as an extended .xyz file, containing one or multiple frames.
        /// Each frame is expected to be in format: <br>
        ///  line 1: N - number of atoms <br>
        ///  line 2: Header containing property map <br>
        ///  lines 3-(N+2): per-atom data
        ///
        /// The property map should be in format: <br>
        ///  prop_name=prop_val : for properties that can be written without whitespaces, <br>
        ///  prop_name="val1 val2" : otherwise <br>
        ///
        /// Allowed types of properties are defined in @ref XYZInfoType,
        /// however, conversion from string will be attempted only if prop_name matches one in main_props,
        /// moreover, an error will be thrown if such conversion is impossible. <br>
        ///
        /// The property name "Properties" is reserved exclusively for array type definition,
        /// and if not found will default to "species:S:1:pos:R:3". <br>
        /// "Properties" defines what array names, types, and in what order they will be expected later lines. <br>
        /// Type support can be found in @ref XYZArrayType, which should be encoded as: <br>
        ///  I:1 - int, <br>
        ///  R:1 - @ref Real, <br>
        ///  R:3 - @ref Vector3, <br>
        ///  S:1 - string without whitespaces if (array name != main_props.species), @ref Species otherwise.
        ///
        ///  PBC is expected to be defined as pbc="T T T";
        ///  pbc="true true true" and pbc="1 1 1" are allowed but will fall back to pbc="T T T" upon writing. <br>
        ///  Lattice can be defined as: "a1 a2 a3 b1 b2 b3 c1 c2 c3", or "a1 b2 c3" if off-diagonals are zero. <br>
        ///  Virial stress tensor must be defined fully "xx xy xz yx yy yz zx zy zz",
        ///  and in addition it is checked that the expected symmetry is maintained up to a numeric error. <br>
        ///
        /// @param in_stream input stream.
        /// @param main_props labels expected for some basic properties.
        /// @return XYZ frames.
        static std::vector<XYZData> read(std::istream& in_stream, const MainXYZPropertyNames& main_props = {});

        XYZData() = default;
        XYZData(const MainXYZPropertyNames& main_props) : main_property_names(main_props) {};
        XYZData(const XYZData& other) = default;
        XYZData(XYZData&& other) = default;
        XYZData& operator=(XYZData&& other) = default;
        XYZData& operator=(const XYZData& other) = default;

        /// Opposite of @ref XYZData::read.
        /// Somewhat more permissive in terms of property names.
        void write(const std::string& filename) const;
        void write(std::ostream& out_stream) const;
        std::string write() const;

        const MainXYZPropertyNames& getMainPropertyNames() const { return main_property_names; }

        const auto& getProperties() const { return info; }
        const auto& getArrays() const { return arrays; }

        std::map<std::string, XYZInfoType>& getPropertiesForEditing();
        std::map<std::string, XYZArrayType>& getArraysForEditing();

    protected:
        MainXYZPropertyNames main_property_names;

        std::map<std::string, XYZInfoType> info;
        std::map<std::string, XYZArrayType> arrays;
    };
}

#endif
