#include "utils/Utils.hpp"

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>
#include <sstream>
#include <fstream>
#include <ranges>
#include <Eigen/Dense>
#include <deque>
#include <unistd.h>

#include "core/neighbours/NeighbourFinder.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {
    std::map<std::string, std::string> parseHeaderLine(const std::string &line) {
        std::map<std::string, std::string> header;

        try {
            size_t pos = 0;
            while (pos < line.size()) {
                if (isspace(line[pos])) {
                    pos++;
                    continue;
                }

                std::string property;
                while (line[pos] != '=') {
                    property += line[pos];
                    pos++;

                    if (pos >= line.size() || isspace(line[pos])) {
                        throw std::runtime_error("'=' not found after " + property);
                    }
                }
                pos++;

                std::string value;
                if (line[pos] == '"') {
                    pos++;
                    while (pos < line.size() && line[pos] != '"') {
                        value += line[pos];
                        pos++;
                    }
                } else {
                    while (pos < line.size() && !isspace(line[pos])) {
                        value += line[pos];
                        pos++;
                    }
                }
                pos++;

                header[property] = value;
            }
        } catch (std::exception& e) {
            JGAP_LOG_AND_THROW("Formatting error {} in : {}", e.what(), line);
        } catch (...) {
            JGAP_LOG_AND_THROW("Formatting error in: {}", line);
        }


        return header;
    }

    bool getLine(std::ifstream &file, std::string &line) {
        if (!getline(file, line)) return false;
        if (!line.empty() && line.back() == '\r') {
            line.pop_back(); // remove Windows carriage return
        }
        return true;
    }

    std::string uniqueStamp() {
        using namespace std::chrono;

        auto now = system_clock::now();
        auto tt = system_clock::to_time_t(now);

        std::tm local_tm{};
#if defined(_WIN32)
        localtime_s(&local_tm, &tt);
#else
        localtime_r(&tt, &local_tm);
#endif

        auto ms = duration_cast<milliseconds>(now.time_since_epoch()) % 1000;

        std::ostringstream oss;
        // Human-readable and filename-friendly: YYYY-MM-DD_HH-MM-SS_mmm-pPID
        oss << std::put_time(&local_tm, "%Y-%m-%d_%H-%M-%S")
            << '_' << std::setw(3) << std::setfill('0') << ms.count();

#if defined(_WIN32)
        DWORD pid = GetCurrentProcessId();
        oss << "-p" << pid;
#else
        oss << "-p" << getpid();
#endif

        return oss.str();
    }

    double factorial(size_t n) {
        double result = 1.0;
        for (int i = 2; i <= n; ++i)
            result *= i;
        return result;
    }

    std::vector<AtomicStructure> readXyz(const std::string& fileName) {

        std::vector<AtomicStructure> result;
        std::ifstream file(fileName);

        if (!file) {
            JGAP_LOG_ERROR( format("Error opening file {}", fileName), true);
        }

        std::string line;
        while (getLine(file, line)) {
            // n_atoms
            size_t n;
            std::istringstream iss(line);
            iss >> n;
            if (!iss.eof()) {
                JGAP_LOG_AND_THROW("Expected single integer, got {}", line);
            }

            // metadata
            getLine(file, line);

            std::map<std::string, std::string> properties = parseHeaderLine(line);

            if (!properties.contains("pbc") || properties["pbc"] != "T T T") {
                JGAP_LOG_AND_THROW("No PBC? : {}", line);
            }
            properties.erase("pbc");

            std::array<Vector3, 3> lattice{};
            if (!properties.contains("Lattice")) {
                JGAP_LOG_AND_THROW("Lattice unspecified in {}", line);
            }
            iss = std::istringstream(properties["Lattice"]);
            iss >> lattice[0].x >> lattice[0].y >> lattice[0].z
                >> lattice[1].x >> lattice[1].y >> lattice[1].z
                >> lattice[2].x >> lattice[2].y >> lattice[2].z;
            properties.erase("Lattice");

            std::optional<double> energy{};
            if (properties.contains("energy")) {
                iss = std::istringstream(properties["energy"]);
                double energyVal;
                iss >> energyVal;
                energy = energyVal;
                properties.erase("energy");
            }

            std::optional<double> energySigmaInverse{};
            if (properties.contains("energy_sigma")) {
                iss = std::istringstream(properties["energy_sigma"]);
                double energySigmaVal;
                iss >> energySigmaVal;
                energySigmaInverse = 1.0 / energySigmaVal;
                properties.erase("energy_sigma");
            }

            std::optional<std::array<Vector3, 3>> virials{};
            if (properties.contains("virials")) {
                properties["virial"] = properties["virials"];
                properties.erase("virials");
            }
            if (properties.contains("virial")) {
                iss = std::istringstream(properties["virial"]);
                std::array<Vector3, 3> virialsVal{};
                iss >> virialsVal[0].x >> virialsVal[0].y >> virialsVal[0].z
                    >> virialsVal[1].x >> virialsVal[1].y >> virialsVal[1].z
                    >> virialsVal[2].x >> virialsVal[2].y >> virialsVal[2].z;
                virials = virialsVal;
                properties.erase("virial");
            }

            std::optional<std::array<Vector3, 3>> virialsSigmasInverse{};
            if (properties.contains("virials_sigma")) {
                iss = std::istringstream(properties["virials_sigma"]);
                double virialsSigmaVal;
                iss >> virialsSigmaVal;
                virialsSigmasInverse = std::array{
                    Vector3{1.0 / virialsSigmaVal, 1.0 / virialsSigmaVal, 1.0 / virialsSigmaVal},
                    Vector3{1.0 / virialsSigmaVal, 1.0 / virialsSigmaVal, 1.0 / virialsSigmaVal},
                    Vector3{1.0 / virialsSigmaVal, 1.0 / virialsSigmaVal, 1.0 / virialsSigmaVal}
                };
                properties.erase("virials_sigma");
            } else if (properties.contains("virials_sigmas")) {
                iss = std::istringstream(properties["virials_sigmas"]);
                std::array<Vector3, 3> virialsSigmasVal{};
                iss >> virialsSigmasVal[0].x >> virialsSigmasVal[0].y >> virialsSigmasVal[0].z
                    >> virialsSigmasVal[1].x >> virialsSigmasVal[1].y >> virialsSigmasVal[1].z
                    >> virialsSigmasVal[2].x >> virialsSigmasVal[2].y >> virialsSigmasVal[2].z;
                virialsSigmasInverse = std::array{
                    Vector3{1.0 / virialsSigmasVal[0].x, 1.0 / virialsSigmasVal[0].y, 1.0 / virialsSigmasVal[0].z},
                    Vector3{1.0 / virialsSigmasVal[1].x, 1.0 / virialsSigmasVal[1].y, 1.0 / virialsSigmasVal[1].z},
                    Vector3{1.0 / virialsSigmasVal[2].x, 1.0 / virialsSigmasVal[2].y, 1.0 / virialsSigmasVal[2].z}
                };
                properties.erase("virials_sigmas");
            }

            std::vector<Vector3> positions(n);
            std::optional<std::vector<Vector3>> forces;
            std::optional<std::vector<Vector3>> forceSigmas;
            std::vector<Species> species(n);
            // todo
            if (line.contains("Properties=species:S:1:pos:R:3:force:R:3:force_sigma:R:3")) {
                forces = std::vector<Vector3>(n);
                forceSigmas = std::vector<Vector3>(n);
                for (size_t i = 0; i < n; i++) {
                    getline(file, line);
                    iss = std::istringstream(line);
                    iss >> species[i];
                    iss >> positions[i].x >> positions[i].y >> positions[i].z;
                    iss >> (*forces)[i].x >> (*forces)[i].y >> (*forces)[i].z;
                    iss >> (*forceSigmas)[i].x >> (*forceSigmas)[i].y >> (*forceSigmas)[i].z;
                }
            } else if (line.contains("Properties=species:S:1:pos:R:3:force:R:3")) {
                forces = std::vector<Vector3>(n);
                for (size_t i = 0; i < n; i++) {
                    getline(file, line);
                    iss = std::istringstream(line);
                    iss >> species[i];
                    iss >> positions[i].x >> positions[i].y >> positions[i].z;
                    iss >> (*forces)[i].x >> (*forces)[i].y >> (*forces)[i].z;
                }
            } else if (line.contains("Properties=species:S:1:pos:R:3")) {
                for (size_t i = 0; i < n; i++) {
                    getline(file, line);
                    iss = std::istringstream(line);
                    iss >> species[i];
                    iss >> positions[i].x >> positions[i].y >> positions[i].z;
                }
            } else {
                JGAP_LOG_AND_THROW("Unknown properties string: {}", line);
            }

            result.push_back(AtomicStructure{
                .properties = properties,
                .lattice_vectors = lattice,
                .positions = positions,
                .species = species,
                .energy = energy,
                .forces = forces,
                .virials = virials,
                .energy_sigma_inverse = energySigmaInverse,
                .force_sigmas_inverse = forceSigmas.transform([](std::vector<Vector3> v) {
                    return v | std::views::transform([](const Vector3& v_i) {
                        return Vector3{1.0 / v_i.x,1.0 / v_i.x,1.0 / v_i.x};
                    }) | std::ranges::to<std::vector>();
                }),
                .virial_sigmas_inverse = virialsSigmasInverse
            });
        }

        return result;
    }

    std::vector<AtomicStructure> readXyz(const std::string &fileName, const double cutoff) {
        auto result = readXyz(fileName);
        NeighbourFinder::findNeighbours(result, cutoff);
        return result;
    }

    void writeXyz(const std::string &fileName, const std::vector<AtomicStructure> &structures) {
        std::ofstream outputStream(fileName);
        writeXyz(outputStream, structures);
        outputStream.close();
    }

    void writeXyz(std::ofstream &outputStream, const std::vector<AtomicStructure> &structures) {
        for (auto& structure: structures) {
            outputStream << structure.size() << std::endl;

            std::string meta;
            meta += "pbc=\"T T T\" ";
            meta += "Lattice=\"";
            meta += std::format(
                "{} {} {} {} {} {} {} {} {}",
                structure.lattice_vectors[0].x, structure.lattice_vectors[0].y, structure.lattice_vectors[0].z,
                structure.lattice_vectors[1].x, structure.lattice_vectors[1].y, structure.lattice_vectors[1].z,
                structure.lattice_vectors[2].x, structure.lattice_vectors[2].y, structure.lattice_vectors[2].z
                );
            meta += "\" ";
            for (auto &[k, v]: structure.properties) {
                if (v.contains(" ") || v.contains("\t")) {
                   meta += k + "=\"" + v + "\" ";
                } else {
                    meta += k + "=" + v + " ";
                }
            }
            if (structure.energy.has_value()) {
                meta += "energy=" + std::to_string(structure.energy.value()) + " ";
            }
            if (structure.virials.has_value()) {
                meta += "virials=\"";
                meta += std::format(
                    "{} {} {} {} {} {} {} {} {}",
                    structure.virials.value()[0].x, structure.virials.value()[0].y, structure.virials.value()[0].z,
                    structure.virials.value()[1].x, structure.virials.value()[1].y, structure.virials.value()[1].z,
                    structure.virials.value()[2].x, structure.virials.value()[2].y, structure.virials.value()[2].z
                );
                meta += "\" ";
            }
            if (structure.forces.has_value()) {
                meta += "Properties=species:S:1:pos:R:3:force:R:3";
            } else {
                meta += "Properties=species:S:1:pos:R:3";
            }

            outputStream << meta << std::endl;

            for (const auto& atom: structure) {
                outputStream << atom.species() << " ";
                outputStream << atom.position().x << " " << atom.position().y << " " << atom.position().z << " ";
                if (structure.forces.has_value()) {
                    outputStream << atom.force().x << " " << atom.force().y << " " << atom.force().z;
                }
                outputStream << std::endl;
            }
        }
        outputStream.flush();
    }

    double rms(const std::vector<double> &x) {
        if (x.empty()) return 0.0;
        double res = 0.0;
        for (double i : x) {
            res += i*i;
        }
        res = sqrt(res / static_cast<double>(x.size()));
        return res;
    }

    double deviation(const std::vector<double> &x) {
        if (x.empty()) return 0.0;
        double mean = 0.0;
        for (double i : x) {
            mean += i;
        }
        mean /= static_cast<double>(x.size());
        double variance = 0.0;
        for (double i : x) {
            variance += (i - mean) * (i - mean);
        }
        variance = variance / static_cast<double>(x.size());
        return sqrt(variance);
    }

    std::vector<std::string> split(const std::string& s, char delimiter) {
        std::vector<std::string> result;
        std::stringstream ss(s);
        std::string token;

        while (getline(ss, token, delimiter)) {
            result.push_back(token);
        }

        return result;
    }

    std::string join(const std::vector<std::string> &strs, char delimiter) {
        std::string res = strs[0];
        for (size_t i = 1; i < strs.size(); i++) {
            res += std::string{delimiter} + strs[i];
        }
        return res;
    }

    std::string withoutExtension(const std::string &s) {

        if (s.find('.') == std::string::npos) {
            return s;
        }

        auto after_split = split(s, '.');
        return join(std::vector(after_split.begin(), after_split.end() - 1), '.');
    }

    std::string matrixToString(const Eigen::MatrixXd& mat) {
        std::stringstream ss;
        for (int i = 0; i < mat.rows(); ++i) {
            for (int j = 0; j < mat.cols(); ++j) {
                ss << mat(i, j);
                if (j != mat.cols() - 1)
                    ss << ", ";
            }
            ss << "\n";
        }
        return ss.str();
    }

    std::string vectorToString(const Eigen::VectorXd &vec) {
        std::stringstream ss;
        for (int i = 0; i < vec.size(); ++i) {
            ss << vec[i];
            if (i != vec.size() - 1) {
                ss << ",";
            }
        }
        return ss.str();
    }

    std::string vectorToString(const std::vector<double> &vec) {
        std::stringstream ss;
        for (int i = 0; i < vec.size(); ++i) {
            ss << vec[i];
            if (i != vec.size() - 1) {
                ss << ",";
            }
        }
        return ss.str();
    }

    std::string vectorToString(const std::vector<size_t> &vec) {
        std::stringstream ss;
        for (int i = 0; i < vec.size(); ++i) {
            ss << vec[i];
            if (i != vec.size() - 1) {
                ss << ",";
            }
        }
        return ss.str();
    }

    std::string vectorToString(const std::vector<std::string> &vec) {
        std::stringstream ss;
        for (int i = 0; i < vec.size(); ++i) {
            ss << vec[i];
            if (i != vec.size() - 1) {
                ss << ",";
            }
        }
        return ss.str();
    }

    static bool nodeIsObject(const DataNode& n) {
        return n.type == DataNode::Type::OBJECT;
    }

    static bool nodeIsArray(const DataNode& n) {
        return n.type == DataNode::Type::ARRAY;
    }

    DataNode& requireFull(DataNode& n,
                          const std::string& key,
                          const char* file,
                          int line,
                          const char* function) {
        if (!nodeIsObject(n)) {
            std::ostringstream oss;
            oss << "Node is not an object when accessing key: \"" << key << "\"\n"
                << "  at " << file << ":" << line << "\n"
                << "  in " << function;
            throw std::domain_error(oss.str());
        }
        auto& m = std::get<std::map<std::string, DataNode>>(n.value);
        auto it = m.find(key);
        if (it == m.end()) {
            std::ostringstream oss;
            oss << "Node key not found: \"" << key << "\"\n"
                << "  at " << file << ":" << line << "\n"
                << "  in " << function;
            throw std::out_of_range(oss.str());
        }
        return it->second;
    }

    const DataNode& requireFull(const DataNode& n,
                                const std::string& key,
                                const char* file,
                                int line,
                                const char* function) {
        if (!nodeIsObject(n)) {
            std::ostringstream oss;
            oss << "DataNode is not an object when accessing key: \"" << key << "\"\n"
                << "  at " << file << ":" << line << "\n"
                << "  in " << function;
            throw std::domain_error(oss.str());
        }
        const auto& m = std::get<std::map<std::string, DataNode>>(n.value);
        auto it = m.find(key);
        if (it == m.end()) {
            std::ostringstream oss;
            oss << "Node key not found: \"" << key << "\"\n"
                << "  at " << file << ":" << line << "\n"
                << "  in " << function;
            throw std::out_of_range(oss.str());
        }
        return it->second;
    }

    DataNode requireArrayFull(DataNode &n,
                              const char* file,
                              int line,
                              const char* function) {
        if (!nodeIsArray(n)) {
            std::ostringstream oss;
            oss << "Node element is not an array\n"
                << "  at " << file << ":" << line << "\n"
                << "  in " << function;
            throw std::domain_error(oss.str());
        }
        return n;
    }

    const DataNode& requireArrayFull(const DataNode &n,
                                     const char* file,
                                     int line,
                                     const char* function) {
        if (!nodeIsArray(n)) {
            std::ostringstream oss;
            oss << "Node element is not an array\n"
                << "  at " << file << ":" << line << "\n"
                << "  in " << function;
            throw std::domain_error(oss.str());
        }
        return n;
    }
}
