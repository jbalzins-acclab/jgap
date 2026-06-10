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

    struct XYZData {
        std::map<std::string, PropertyValue> properties;
        std::map<std::string, ArrayValue> arrays;

        XYZData() = default;
        XYZData(const XYZData& other) = default;

        XYZData& operator=(const XYZData& other) {
            if (this != &other) {
                properties = other.properties;
                arrays = other.arrays;
            }
            return *this;
        }

        static bool getLine(std::ifstream &file, std::string &line);

        static std::vector<XYZData> read(const std::string &filename) {
            std::ifstream in(filename);

            if (!in.is_open()) {
                JGAP_LOG_AND_THROW("Failed to open file: {}", filename);
            }

            return read(in);
        }
        static std::vector<XYZData> read(std::ifstream &in_stream);

        void write(const std::string &filename) const {
            std::ofstream out(filename);
            write(out);
        }
        void write(std::ofstream &out_stream) const;
    };
}

#endif