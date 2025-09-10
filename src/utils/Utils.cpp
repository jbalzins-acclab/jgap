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

#include "core/neighbours/NeighbourFinder.hpp"
#include "io/log/CurrentLogger.hpp"

using namespace std;

namespace jgap {
    map<string, string> parseHeaderLine(const string &line) {
        map<string, string> header;

        try {
            size_t pos = 0;
            while (pos < line.size()) {
                if (isspace(line[pos])) {
                    pos++;
                    continue;
                }

                string property = "";
                while (line[pos] != '=') {
                    property += line[pos];
                    pos++;

                    if (pos >= line.size() || isspace(line[pos])) {
                        throw runtime_error("'=' not found after " + property);
                    }
                }
                pos++;

                string value = "";
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
        } catch (exception& e) {
            CurrentLogger::get()->logAndThrow("Formatting error {} in : {}", e.what(), line);
        } catch (...) {
            CurrentLogger::get()->logAndThrow("Formatting error in: {}", line);
        }


        return header;
    }

    bool getLine(ifstream &file, string &line) {
        if (!getline(file, line)) return false;
        if (!line.empty() && line.back() == '\r') {
            line.pop_back(); // remove Windows carriage return
        }
        return true;
    }

    vector<AtomicStructure> readXyz(const string& fileName) {

        vector<AtomicStructure> result;
        ifstream file(fileName);

        if (!file) {
            CurrentLogger::get()->error( format("Error opening file {}", fileName), true);
        }

        string line;
        while (getLine(file, line)) {
            // n_atoms
            size_t n;
            istringstream iss(line);
            iss >> n;
            if (!iss.eof()) {
                CurrentLogger::get()->logAndThrow("Expected single integer, got {}", line);
            }

            // metadata
            getLine(file, line);

            map<string, string> properties = parseHeaderLine(line);

            if (!properties.contains("pbc") || properties["pbc"] != "T T T") {
                CurrentLogger::get()->logAndThrow("No PBC? : {}", line);
            }
            properties.erase("pbc");

            array<Vector3, 3> lattice{};
            if (!properties.contains("Lattice")) {
                CurrentLogger::get()->logAndThrow("Lattice unspecified in {}", line);
            }
            iss = istringstream(properties["Lattice"]);
            iss >> lattice[0].x >> lattice[0].y >> lattice[0].z
                >> lattice[1].x >> lattice[1].y >> lattice[1].z
                >> lattice[2].x >> lattice[2].y >> lattice[2].z;
            properties.erase("Lattice");

            optional<double> energy{};
            if (properties.contains("energy")) {
                iss = istringstream(properties["energy"]);
                double energyVal;
                iss >> energyVal;
                energy = energyVal;
                properties.erase("energy");
            }

            optional<double> energySigmaInverse{};
            if (properties.contains("energy_sigma")) {
                iss = istringstream(properties["energy_sigma"]);
                double energySigmaVal;
                iss >> energySigmaVal;
                energySigmaInverse = 1.0 / energySigmaVal;
                properties.erase("energy_sigma");
            }

            optional<array<Vector3, 3>> virials{};
            if (properties.contains("virials")) {
                properties["virial"] = properties["virials"];
                properties.erase("virials");
            }
            if (properties.contains("virial")) {
                iss = istringstream(properties["virial"]);
                array<Vector3, 3> virialsVal{};
                iss >> virialsVal[0].x >> virialsVal[0].y >> virialsVal[0].z
                    >> virialsVal[1].x >> virialsVal[1].y >> virialsVal[1].z
                    >> virialsVal[2].x >> virialsVal[2].y >> virialsVal[2].z;
                virials = virialsVal;
                properties.erase("virial");
            }

            optional<array<Vector3, 3>> virialsSigmasInverse{};
            if (properties.contains("virials_sigma")) {
                iss = istringstream(properties["virials_sigma"]);
                double virialsSigmaVal;
                iss >> virialsSigmaVal;
                virialsSigmasInverse = array{
                    Vector3{1.0 / virialsSigmaVal, 1.0 / virialsSigmaVal, 1.0 / virialsSigmaVal},
                    Vector3{1.0 / virialsSigmaVal, 1.0 / virialsSigmaVal, 1.0 / virialsSigmaVal},
                    Vector3{1.0 / virialsSigmaVal, 1.0 / virialsSigmaVal, 1.0 / virialsSigmaVal}
                };
                properties.erase("virials_sigma");
            } else if (properties.contains("virials_sigmas")) {
                iss = istringstream(properties["virials_sigmas"]);
                array<Vector3, 3> virialsSigmasVal{};
                iss >> virialsSigmasVal[0].x >> virialsSigmasVal[0].y >> virialsSigmasVal[0].z
                    >> virialsSigmasVal[1].x >> virialsSigmasVal[1].y >> virialsSigmasVal[1].z
                    >> virialsSigmasVal[2].x >> virialsSigmasVal[2].y >> virialsSigmasVal[2].z;
                virialsSigmasInverse = array{
                    Vector3{1.0 / virialsSigmasVal[0].x, 1.0 / virialsSigmasVal[0].y, 1.0 / virialsSigmasVal[0].z},
                    Vector3{1.0 / virialsSigmasVal[1].x, 1.0 / virialsSigmasVal[1].y, 1.0 / virialsSigmasVal[1].z},
                    Vector3{1.0 / virialsSigmasVal[2].x, 1.0 / virialsSigmasVal[2].y, 1.0 / virialsSigmasVal[2].z}
                };
                properties.erase("virials_sigmas");
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
                CurrentLogger::get()->logAndThrow("Unknown properties string: {}", line);
            }

            result.push_back(AtomicStructure{
                .properties = properties,
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

    vector<AtomicStructure> readXyz(const string &fileName, const double cutoff) {
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
            for (auto &[k, v]: structure.properties) {
                if (v.contains(" ") || v.contains("\t")) {
                   meta += k + "=\"" + v + "\"";
                } else {
                    meta += k + "=" + v;
                }
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
