#ifndef JGAP_XYZDATA_HPP
#define JGAP_XYZDATA_HPP
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <variant>

#include "core/Real.hpp"
#include "core/atomic/geometry/Vector3.hpp"
#include "core/atomic/energy/Virials.hpp"
#include "core/atomic/geometry/Lattice.hpp"
#include "core/atomic/species/Species.hpp"

namespace jgap {
    using PropertyValue = std::variant<std::string, int, Real, Vector3, Virials, Lattice>;
    using ArrayValue = std::variant<
        std::vector<int>, std::vector<Real>, std::vector<Vector3>, std::vector<std::string>, std::vector<Species>
    >;

    struct XYZData {
        std::map<std::string, PropertyValue> properties;
        std::map<std::string, ArrayValue> arrays;

        static XYZData read(const std::string &filename) {
            std::ifstream in(filename);
            return read(in);
        }
        static XYZData read(std::ifstream &in_stream);

        void write(const std::string &filename) const {
            std::ofstream out(filename);
            write(out);
        }
        void write(std::ofstream &out_stream) const;
    };
}

#endif
