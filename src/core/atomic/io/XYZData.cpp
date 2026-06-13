#include "XYZData.hpp"
#include "utils/Utils.hpp"
#include "core/atomic/geometry/Lattice.hpp"
#include "core/atomic/energy/Virials.hpp"
#include "core/atomic/species/Species.hpp"
#include "io/log/CurrentLogger.hpp"
#include <sstream>
#include <iomanip>
#include <variant>

namespace jgap {

    std::vector<XYZData> XYZData::read(const std::string &filename) {
        std::ifstream in(filename);

        if (!in.is_open()) {
            JGAP_LOG_AND_THROW("Failed to open file: {}", filename);
        }

        return read(in);
    }

    std::vector<XYZData> XYZData::read(std::ifstream &in_stream) {
        std::vector<XYZData> frames;

        while (true) {
            std::string line;
            if (!getLine(in_stream, line)) break;

            size_t n_atoms = 0;
            try {
                n_atoms = std::stoul(line);
            } catch (...) {
                break;
            }

            XYZData data;
            if (!getLine(in_stream, line)) break;

            auto raw_properties = parseHeaderLine(line);

            for (auto const& [k, v] : raw_properties) {
                if (k == "Properties") continue;

                if (k == "pbc") {
                    std::istringstream iss(v);
                    std::string token;
                    std::array pbc_vals = {false, false, false};
                    int idx = 0;
                    while (iss >> token && idx < 3) {
                        pbc_vals[idx++] = (token == "T" || token == "true" || token == "1");
                    }
                    data.properties[k] = pbc_vals;
                    continue;
                }

                std::istringstream iss(v);
                std::vector<Real> vals;
                Real val;
                while (iss >> val) {
                    vals.push_back(val);
                }

                if (k == "Lattice" || k == "lattice") {
                    if (vals.size() == 9) {
                        data.properties[k] = Lattice{
                            Vector3{vals[0], vals[1], vals[2]},
                            Vector3{vals[3], vals[4], vals[5]},
                            Vector3{vals[6], vals[7], vals[8]}
                        };
                    } else if (vals.size() == 3) {
                        data.properties[k] = Lattice{
                            Vector3{vals[0], 0, 0},
                            Vector3{0, vals[1], 0},
                            Vector3{0, 0, vals[2]}
                        };
                    } else {
                        data.properties[k] = v;
                    }
                    continue;
                }

                if (vals.size() == 1) {
                    if (v.find('.') != std::string::npos) {
                        data.properties[k] = (Real)vals[0];
                    } else {
                        try {
                            data.properties[k] = std::stoi(v);
                        } catch (...) {
                            data.properties[k] = (Real)vals[0];
                        }
                    }
                } else if (vals.size() == 3) {
                    data.properties[k] = Vector3{vals[0], vals[1], vals[2]};
                } else if (vals.size() == 9) {
                    if (std::abs(vals[1] - vals[3]) < 1e-9 && std::abs(vals[2] - vals[6]) < 1e-9 && std::abs(vals[5] - vals[7]) < 1e-9) {
                         data.properties[k] = Virials{vals[0], vals[1], vals[2], vals[4], vals[5], vals[8]};
                    } else {
                         data.properties[k] = v;
                    }
                } else {
                    data.properties[k] = v;
                }
            }

            struct PropInfo {
                std::string name;
                char type;
                int count;
            };
            std::vector<PropInfo> prop_infos;

            if (raw_properties.contains("Properties")) {
                std::string props_str = raw_properties["Properties"];
                auto tokens = split(props_str, ':');
                for (size_t i = 0; i < tokens.size(); i += 3) {
                    prop_infos.push_back({tokens[i], tokens[i+1][0], std::stoi(tokens[i+2])});
                }
            } else {
                prop_infos.push_back({"species", 'S', 1});
                prop_infos.push_back({"pos", 'R', 3});
            }

            for (size_t i = 0; i < n_atoms; ++i) {
                if (!getLine(in_stream, line)) {
                    JGAP_LOG_AND_THROW("Unexpected end of file at atom " + std::to_string(i));
                }
                std::istringstream iss(line);
                for (const auto& info : prop_infos) {
                    if (info.type == 'R') {
                        if (info.count == 3) {
                            Vector3 v;
                            if (!(iss >> v.x >> v.y >> v.z)) JGAP_LOG_AND_THROW("Failed to read Vector3 for property " + info.name + " at atom " + std::to_string(i));
                            if (!data.arrays.contains(info.name)) data.arrays[info.name] = std::vector<Vector3>();
                            std::get<std::vector<Vector3>>(data.arrays[info.name]).push_back(v);
                        } else {
                            if (!data.arrays.contains(info.name)) data.arrays[info.name] = std::vector<Real>();
                            for (int c = 0; c < info.count; ++c) {
                                Real r;
                                if (!(iss >> r)) JGAP_LOG_AND_THROW("Failed to read Real for property " + info.name + " at atom " + std::to_string(i));
                                std::get<std::vector<Real>>(data.arrays[info.name]).push_back(r);
                            }
                        }
                    } else if (info.type == 'I') {
                        if (!data.arrays.contains(info.name)) data.arrays[info.name] = std::vector<int>();
                        for (int c = 0; c < info.count; ++c) {
                            int val;
                            if (!(iss >> val)) JGAP_LOG_AND_THROW("Failed to read int for property " + info.name + " at atom " + std::to_string(i));
                            std::get<std::vector<int>>(data.arrays[info.name]).push_back(val);
                        }
                    } else if (info.type == 'S') {
                        if (info.name == "species") {
                            if (!data.arrays.contains(info.name)) data.arrays[info.name] = std::vector<Species>();
                            for (int c = 0; c < info.count; ++c) {
                                std::string s;
                                if (!(iss >> s)) JGAP_LOG_AND_THROW("Failed to read string for property " + info.name + " at atom " + std::to_string(i));
                                std::get<std::vector<Species>>(data.arrays[info.name]).push_back(Species(s));
                            }
                        } else {
                            if (!data.arrays.contains(info.name)) data.arrays[info.name] = std::vector<std::string>();
                            for (int c = 0; c < info.count; ++c) {
                                std::string s;
                                if (!(iss >> s)) JGAP_LOG_AND_THROW("Failed to read string for property " + info.name + " at atom " + std::to_string(i));
                                std::get<std::vector<std::string>>(data.arrays[info.name]).push_back(s);
                            }
                        }
                    }
                }
            }

            frames.push_back(std::move(data));
        }

        return frames;
    }

    XYZData& XYZData::operator=(const XYZData &other) {
        if (this != &other) {
            properties = other.properties;
            arrays = other.arrays;
        }
        return *this;
    }

    void XYZData::write(const std::string &filename) const {
        std::ofstream out(filename);
        write(out);
    }

    void XYZData::write(std::ofstream &out_stream) const {
        size_t n_atoms = 0;
        for (auto const& [name, val] : arrays) {
             size_t current_size = 0;
             std::visit([&current_size](auto&& arg) { current_size = arg.size(); }, val);
             if (n_atoms == 0) n_atoms = current_size;
             else if (n_atoms != current_size) {
                 JGAP_LOG_AND_THROW("Array sizes mismatch in XYZData::write: " + name + " has size " + std::to_string(current_size) + " but expected " + std::to_string(n_atoms));
             }
        }

        out_stream << n_atoms << "\n";

        std::vector<std::string> meta_parts;
        for (auto const& [k, v] : properties) {
            std::string key = k;
            if (std::holds_alternative<Lattice>(v)) key = "Lattice";

            std::string part = key + "=";
            std::string val_str;
            std::visit([&val_str](auto&& arg) {
                using T = std::decay_t<decltype(arg)>;
                if constexpr (std::is_same_v<T, std::string>) val_str = arg;
                else if constexpr (std::is_same_v<T, int>) val_str = std::to_string(arg);
                else if constexpr (std::is_same_v<T, Real>) val_str = std::to_string(arg);
                else if constexpr (std::is_same_v<T, Vector3>) val_str = std::to_string(arg.x) + " " + std::to_string(arg.y) + " " + std::to_string(arg.z);
                else if constexpr (std::is_same_v<T, Virials>) {
                    std::ostringstream oss;
                    oss << arg.xx << " " << arg.xy << " " << arg.xz << " " << arg.xy << " " << arg.yy << " " << arg.yz << " " << arg.xz << " " << arg.yz << " " << arg.zz;
                    val_str = oss.str();
                } else if constexpr (std::is_same_v<T, Lattice>) {
                    std::ostringstream oss;
                    oss << arg.a.x << " " << arg.a.y << " " << arg.a.z << " " << arg.b.x << " " << arg.b.y << " " << arg.b.z << " " << arg.c.x << " " << arg.c.y << " " << arg.c.z;
                    val_str = oss.str();
                } else if constexpr (std::is_same_v<T, std::array<bool, 3>>) {
                    std::ostringstream oss;
                    oss << (arg[0] ? "T" : "F") << " " << (arg[1] ? "T" : "F") << " " << (arg[2] ? "T" : "F");
                    val_str = oss.str();
                }
            }, v);

            if (val_str.find(' ') != std::string::npos || val_str.empty()) part += "\"" + val_str + "\"";
            else part += val_str;
            meta_parts.push_back(part);
        }

        std::vector<std::string> prop_tokens;
        if (arrays.contains("species")) prop_tokens.push_back("species:S:1");
        if (arrays.contains("pos")) prop_tokens.push_back("pos:R:3");

        for (const auto& [name, val] : arrays) {
            if (name == "species" || name == "pos") continue;
            std::visit([&prop_tokens, &name](auto&& arg) {
                using T = std::decay_t<decltype(arg)>;
                if constexpr (std::is_same_v<T, std::vector<int>>) prop_tokens.push_back(name + ":I:1");
                else if constexpr (std::is_same_v<T, std::vector<Real>>) prop_tokens.push_back(name + ":R:1");
                else if constexpr (std::is_same_v<T, std::vector<Vector3>>) prop_tokens.push_back(name + ":R:3");
                else if constexpr (std::is_same_v<T, std::vector<std::string>>) prop_tokens.push_back(name + ":S:1");
                else if constexpr (std::is_same_v<T, std::vector<Species>>) prop_tokens.push_back(name + ":S:1");
            }, val);
        }

        meta_parts.push_back("Properties=" + join(prop_tokens, ':'));
        out_stream << join(meta_parts, ' ') << "\n";

        for (size_t i = 0; i < n_atoms; ++i) {
            std::vector<std::string> line_tokens;
            if (arrays.contains("species")) {
                auto const& spec_val = arrays.at("species");
                if (std::holds_alternative<std::vector<Species>>(spec_val)) {
                    line_tokens.push_back(std::get<std::vector<Species>>(spec_val)[i].symbol());
                } else {
                    line_tokens.push_back(std::get<std::vector<std::string>>(spec_val)[i]);
                }
            }
            if (arrays.contains("pos")) {
                const Vector3& pos = std::get<std::vector<Vector3>>(arrays.at("pos"))[i];
                line_tokens.push_back(std::to_string(pos.x));
                line_tokens.push_back(std::to_string(pos.y));
                line_tokens.push_back(std::to_string(pos.z));
            }
            for (const auto& [name, val] : arrays) {
                if (name == "species" || name == "pos") continue;
                std::visit([&line_tokens, i](auto&& arg) {
                    using T = std::decay_t<decltype(arg)>;
                    if constexpr (std::is_same_v<T, std::vector<Vector3>>) {
                        line_tokens.push_back(std::to_string(arg[i].x));
                        line_tokens.push_back(std::to_string(arg[i].y));
                        line_tokens.push_back(std::to_string(arg[i].z));
                    } else if constexpr (std::is_same_v<T, std::vector<std::string>>) {
                        line_tokens.push_back(arg[i]);
                    } else if constexpr (std::is_same_v<T, std::vector<Species>>) {
                        line_tokens.push_back(static_cast<Species>(arg[i]).symbol());
                    } else {
                        line_tokens.push_back(std::to_string(arg[i]));
                    }
                }, val);
            }
            out_stream << join(line_tokens, ' ') << "\n";
        }
    }
}
