#ifndef JGAP_XYZDATA_HPP
#define JGAP_XYZDATA_HPP
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <variant>
#include <array>

#include "core/Real.hpp"
#include "core/atomic/geometry/Vector3.hpp"
#include "core/atomic/energy/Virials.hpp"
#include "core/atomic/geometry/Lattice.hpp"
#include "core/atomic/species/Species.hpp"

namespace jgap {
    using PropertyValue = std::variant<
        std::string,
        int,
        Real,
        Vector3,
        Virials,
        Lattice,
        std::array<bool, 3>
    >;
    using ArrayValue = std::variant<
        std::vector<int>,
        std::vector<Real>,
        std::vector<Vector3>,
        std::vector<std::string>,
        std::vector<Species>
    >;

    class XYZData {
    public:
        static std::vector<XYZData> read(const std::string &filename);
        static std::vector<XYZData> read(std::ifstream &in_stream);

        std::map<std::string, PropertyValue> properties;
        std::map<std::string, ArrayValue> arrays;

        XYZData() = default;
        XYZData(const XYZData& other) = default;

        XYZData& operator=(const XYZData& other);

        void write(const std::string &filename) const;

        void write(std::ofstream &out_stream) const;
    };
}

#endif