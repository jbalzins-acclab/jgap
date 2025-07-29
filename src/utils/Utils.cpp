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

#include "core/neighbours/NeighbourFinder.hpp"
#include "io/log/CurrentLogger.hpp"

using namespace std;

namespace jgap {

    vector<AtomicStructure> readXyz(const string& fileName) {

        vector<AtomicStructure> result;
        ifstream file(fileName);

        if (!file) {
            CurrentLogger::get()->error( format("Error opening file {}", fileName), true);
        }

        string line;
        istringstream iss;
        while (getline(file, line)) {
            // n_atoms
            size_t n;
            iss = istringstream(line);
            iss >> n;
            if (!iss.eof()) {
                CurrentLogger::get()->error( format("Expected single integer, got {}", line), true);
            }

            // metadata
            getline(file, line);
            line += " "; // simplify end of line parsing

            if (!line.contains("pbc=\"T T T\"")) {
                CurrentLogger::get()->error(format("No PBC? : {}", line), true);
            }

            array<Vector3, 3> lattice{};
            if (size_t latticeStartIdx = line.find("Lattice=\""); latticeStartIdx != string::npos) {
                latticeStartIdx += string("Lattice=\"").size();

                size_t latticeEndIdx = line.find('\"', latticeStartIdx);
                if (latticeEndIdx == string::npos) {
                    CurrentLogger::get()->error(format("Lattice unspecified in {}", line), true);
                }

                string latticeStr = line.substr(latticeStartIdx, latticeEndIdx - latticeStartIdx);
                iss = istringstream(latticeStr);
                iss >> lattice[0].x >> lattice[0].y >> lattice[0].z
                    >> lattice[1].x >> lattice[1].y >> lattice[1].z
                    >> lattice[2].x >> lattice[2].y >> lattice[2].z;
            }
            else {
                CurrentLogger::get()->error(format("Lattice unspecified in {}", line), true);
            }

            optional<string> configType{};
            if (size_t configTypeStartIdx = line.find("config_type="); configTypeStartIdx != string::npos) {
                configTypeStartIdx += string("config_type=").size();
                size_t configTypeEndIdx = line.find(' ', configTypeStartIdx);
                if (configTypeEndIdx == string::npos) {
                    CurrentLogger::get()->error(format("Config type formatting error in: {}", line), true);
                }

                configType = line.substr(configTypeStartIdx, configTypeEndIdx - configTypeStartIdx);
            }

            optional<double> energy{};
            if (size_t energyStartIdx = line.find("energy="); energyStartIdx != string::npos) {
                energyStartIdx += string("energy=").size();
                size_t energyEndIdx = line.find(' ', energyStartIdx);
                if (energyEndIdx == string::npos) {
                    CurrentLogger::get()->error(format("Energy formatting error in: {}", line), true);
                }

                string energyStr = line.substr(energyStartIdx, energyEndIdx - energyStartIdx);
                iss = istringstream(energyStr);
                double energyVal;
                iss >> energyVal;
                energy = energyVal;
            }

            optional<double> energySigmaInverse{};
            if (size_t energySigmaStartIdx = line.find("energy_sigma="); energySigmaStartIdx != string::npos) {
                energySigmaStartIdx += string("energy_sigma=").size();
                size_t energySigmaEndIdx = line.find(' ', energySigmaStartIdx);
                if (energySigmaEndIdx == string::npos) {
                    CurrentLogger::get()->error(format("energy_sigma formatting error in: {}", line), true);
                }

                string energySigmaStr = line.substr(energySigmaStartIdx, energySigmaEndIdx - energySigmaStartIdx);
                iss = istringstream(energySigmaStr);
                double energySigmaVal;
                iss >> energySigmaVal;
                energySigmaInverse = 1.0 / energySigmaVal;
            }

            optional<array<Vector3, 3>> virials{};
            size_t virialsStartIdx = line.find("virial=\"");
            if (virialsStartIdx == string::npos) {
                virialsStartIdx = line.find("virials=\"");
                virialsStartIdx += string("virials=\"").size();
            } else {
                virialsStartIdx += string("virial=\"").size();
            }
            if (virialsStartIdx != string::npos) {

                size_t virialsEndIdx = line.find('\"', virialsStartIdx);
                if (virialsEndIdx == string::npos) {
                    CurrentLogger::get()->error(format("Virials parsing error in: {}", line), true);
                }

                string virialsStr = line.substr(virialsStartIdx, virialsEndIdx - virialsStartIdx);
                iss = istringstream(virialsStr);
                array<Vector3, 3> virialsVal{};
                iss >> virialsVal[0].x >> virialsVal[0].y >> virialsVal[0].z
                    >> virialsVal[1].x >> virialsVal[1].y >> virialsVal[1].z
                    >> virialsVal[2].x >> virialsVal[2].y >> virialsVal[2].z;
                virials = virialsVal;
            }

            optional<array<Vector3, 3>> virialsSigmasInverse{};
            if (size_t virialsSigmaStartIdx = line.find("virials_sigma="); virialsSigmaStartIdx != string::npos) {
                virialsSigmaStartIdx += string("virials_sigma=").size();
                size_t virialsSigmaEndIdx = line.find(' ', virialsSigmaStartIdx);
                if (virialsSigmaEndIdx == string::npos) {
                    CurrentLogger::get()->error(format("virials_sigma formatting error in: {}", line), true);
                }

                string virialsSigmaStr = line.substr(virialsSigmaStartIdx, virialsSigmaEndIdx - virialsSigmaStartIdx);
                iss = istringstream(virialsSigmaStr);
                double virialsSigmaVal;
                iss >> virialsSigmaVal;
                virialsSigmasInverse = array{
                    Vector3{1.0 / virialsSigmaVal, 1.0 / virialsSigmaVal, 1.0 / virialsSigmaVal},
                    Vector3{1.0 / virialsSigmaVal, 1.0 / virialsSigmaVal, 1.0 / virialsSigmaVal},
                    Vector3{1.0 / virialsSigmaVal, 1.0 / virialsSigmaVal, 1.0 / virialsSigmaVal}
                };
            } else if (size_t virialsSigmasStartIdx = line.find("virials_sigmas=\"");
                       virialsSigmasStartIdx != string::npos) {
                virialsSigmasStartIdx += string("virials_sigmas=").size();
                size_t virialsSigmasEndIdx = line.find('\"', virialsSigmasStartIdx);
                if (virialsSigmasEndIdx == string::npos) {
                    CurrentLogger::get()->error(format("Virials sigmas parsing error in: {}", line), true);
                }

                string virialsSigmasStr = line.substr(virialsSigmasStartIdx, virialsSigmasEndIdx - virialsSigmasStartIdx);
                iss = istringstream(virialsSigmasStr);
                array<Vector3, 3> virialsSigmasVal{};
                iss >> virialsSigmasVal[0].x >> virialsSigmasVal[0].y >> virialsSigmasVal[0].z
                    >> virialsSigmasVal[1].x >> virialsSigmasVal[1].y >> virialsSigmasVal[1].z
                    >> virialsSigmasVal[2].x >> virialsSigmasVal[2].y >> virialsSigmasVal[2].z;
                virialsSigmasInverse = array{
                    Vector3{1.0 / virialsSigmasVal[0].x, 1.0 / virialsSigmasVal[0].y, 1.0 / virialsSigmasVal[0].z},
                    Vector3{1.0 / virialsSigmasVal[1].x, 1.0 / virialsSigmasVal[1].y, 1.0 / virialsSigmasVal[1].z},
                    Vector3{1.0 / virialsSigmasVal[2].x, 1.0 / virialsSigmasVal[2].y, 1.0 / virialsSigmasVal[2].z}
                };
            }

            vector<Vector3> positions(n);
            optional<vector<Vector3>> forces;
            optional<vector<Vector3>> forceSigmas;
            vector<Species> species(n);
            // todo
            if (line.contains("Properties=species:S:1:pos:R:3:force:R:3:force_sigma:R:3")) {
                forces = vector<Vector3>(n);
                forceSigmas = vector<Vector3>(n);
                for (size_t i = 0; i < n; i++) {
                    getline(file, line);
                    iss = istringstream(line);
                    iss >> species[i];
                    iss >> positions[i].x >> positions[i].y >> positions[i].z;
                    iss >> (*forces)[i].x >> (*forces)[i].y >> (*forces)[i].z;
                    iss >> (*forceSigmas)[i].x >> (*forceSigmas)[i].y >> (*forceSigmas)[i].z;
                }
            } else if (line.contains("Properties=species:S:1:pos:R:3:force:R:3")) {
                forces = vector<Vector3>(n);
                for (size_t i = 0; i < n; i++) {
                    getline(file, line);
                    iss = istringstream(line);
                    iss >> species[i];
                    iss >> positions[i].x >> positions[i].y >> positions[i].z;
                    iss >> (*forces)[i].x >> (*forces)[i].y >> (*forces)[i].z;
                }
            } else if (line.contains("Properties=species:S:1:pos:R:3")) {
                for (size_t i = 0; i < n; i++) {
                    getline(file, line);
                    iss = istringstream(line);
                    iss >> species[i];
                    iss >> positions[i].x >> positions[i].y >> positions[i].z;
                }
            } else {
                CurrentLogger::get()->error(format("Unknown properties string: {}", line), true);
            }

            result.push_back(AtomicStructure{
                .configType = configType,
                .lattice = lattice,
                .positions = positions,
                .species = species,
                .energy = energy,
                .forces = forces,
                .virials = virials,
                .energySigmaInverse = energySigmaInverse,
                .forceSigmasInverse = forceSigmas.transform([](vector<Vector3> v) {
                    return v | std::views::transform([](const Vector3& v_i) {
                        return Vector3{1.0 / v_i.x,1.0 / v_i.x,1.0 / v_i.x};
                    }) | std::ranges::to<std::vector>();
                }),
                .virialSigmasInverse = virialsSigmasInverse
            });
        }

        return result;
    }

    vector<AtomicStructure> readXyz(const string &fileName, double cutoff) {
        auto result = readXyz(fileName);
        NeighbourFinder::findNeighbours(result, cutoff);
        return result;
    }

    void writeXyz(const string &fileName, const vector<AtomicStructure> &structures) {
        ofstream file(fileName);
        for (auto& structure: structures) {
            file << structure.size() << endl;

            string meta = "";
            meta += "pbc=\"T T T\" ";
            meta += "Lattice=\"";
            meta += format(
                "{} {} {} {} {} {} {} {} {}",
                structure.lattice[0].x, structure.lattice[0].y, structure.lattice[0].z,
                structure.lattice[1].x, structure.lattice[1].y, structure.lattice[1].z,
                structure.lattice[2].x, structure.lattice[2].y, structure.lattice[2].z
                );
            meta += "\" ";
            if (structure.configType.has_value()) {
                meta += "config_type=" + structure.configType.value() + " ";
            }
            if (structure.energy.has_value()) {
                meta += "energy=" + to_string(structure.energy.value()) + " ";
            }
            if (structure.virials.has_value()) {
                meta += "virials=\"";
                meta += format(
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

            file << meta << endl;

            for (const auto& atom: structure) {
                file << atom.species() << " ";
                file << atom.position().x << " " << atom.position().y << " " << atom.position().z << " ";
                if (structure.forces.has_value()) {
                    file << atom.force().x << " " << atom.force().y << " " << atom.force().z;
                }
                file << endl;
            }
        }
        file.flush();
        file.close();
    }

    vector<string> split(const string& s, char delimiter) {
        vector<string> result;
        stringstream ss(s);
        string token;

        while (getline(ss, token, delimiter)) {
            result.push_back(token);
        }

        return result;
    }

    void saveArray(const vector<double> &data, const string &filename) {
        ofstream out(filename, ios::binary);
        out.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(double));
    }

    vector<double> loadArray(const string &filename) {
        ifstream in(filename, ios::binary);
        in.seekg(0, ios::end);
        streamsize size = in.tellg();
        in.seekg(0, ios::beg);

        vector<double> data(size / sizeof(double));
        in.read(reinterpret_cast<char*>(data.data()), size);
        return data;
    }

    string matrixToString(const Eigen::MatrixXd& mat) {
        stringstream ss;
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

    string vectorToString(const Eigen::VectorXd &vec) {
        stringstream ss;
        for (int i = 0; i < vec.size(); ++i) {
            ss << vec[i];
            if (i != vec.size() - 1) {
                ss << ",";
            }
        }
        return ss.str();
    }

    string vectorToString(const vector<double> &vec) {
        stringstream ss;
        for (int i = 0; i < vec.size(); ++i) {
            ss << vec[i];
            if (i != vec.size() - 1) {
                ss << ",";
            }
        }
        return ss.str();
    }

    string vectorToString(const vector<size_t> &vec) {
        stringstream ss;
        for (int i = 0; i < vec.size(); ++i) {
            ss << vec[i];
            if (i != vec.size() - 1) {
                ss << ",";
            }
        }
        return ss.str();
    }
}
