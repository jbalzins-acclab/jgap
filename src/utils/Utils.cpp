#include "utils/Utils.hpp"

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>
#include <string_view>
#include <sstream>
#include <Python.h>
#include <sstream>
#include <fstream>
#include <Eigen/Dense>

using namespace std;

namespace jgap {
    string runBash(const string_view command) {
        static array<char, 1000> buffer{};
        string result;
        const unique_ptr<FILE, decltype(&pclose)> pipe(popen(command.data(), "r"), pclose);
        if (!pipe) {
            throw runtime_error("popen() failed!");
        }
        while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
            result += buffer.data();
        }
        return result;
    }

    string executePythonScript(string_view pythonCode) {

        if (!Py_IsInitialized()) {
            Py_Initialize();
        }

        // Redirect stdout to capture output
        PyRun_SimpleString("import sys, io\nsys.stdout = io.StringIO()");

        // Get __main__ module and its namespace
        PyObject* pModule = PyImport_AddModule("__main__");
        if (!pModule) throw runtime_error("Failed to get __main__ module");

        PyObject* pDict = PyModule_GetDict(pModule);

        // Execute user code
        PyObject* pResult = PyRun_String(pythonCode.data(), Py_file_input, pDict, pDict);
        if (!pResult) {
            PyErr_Print();
            throw runtime_error("Python script execution failed");
        }
        Py_DECREF(pResult);

        // Fetch sys.stdout.getvalue()
        PyObject* sysModule = PyImport_ImportModule("sys");
        if (!sysModule) throw runtime_error("Failed to import sys");

        PyObject* stdoutObj = PyObject_GetAttrString(sysModule, "stdout");
        PyObject* getvalueFunc = PyObject_GetAttrString(stdoutObj, "getvalue");

        if (!stdoutObj || !getvalueFunc || !PyCallable_Check(getvalueFunc)) {
            Py_XDECREF(sysModule);
            Py_XDECREF(stdoutObj);
            Py_XDECREF(getvalueFunc);
            throw runtime_error("Failed to access sys.stdout.getvalue");
        }

        PyObject* outputStr = PyObject_CallObject(getvalueFunc, nullptr);
        if (!outputStr) {
            PyErr_Print();
            throw runtime_error("Failed to call sys.stdout.getvalue()");
        }

        const char* output = PyUnicode_AsUTF8(outputStr);
        string result = output ? output : "";

        // Cleanup
        Py_DECREF(outputStr);
        Py_DECREF(getvalueFunc);
        Py_DECREF(stdoutObj);
        Py_DECREF(sysModule);

        return result;
    }

    string xyzToJgapInput(const string_view xyzFileName) {
        string pythonScript = R"(
import ase.io, os
boxes = ase.io.read('__filename__', index=':')

print(len(boxes))
for box in boxes:
    # Mandatory:
    for i in range(3):
        for j in range(3):
            print(box.get_cell(complete=True)[i][j])
    print(len(box))
    for atom in box:
        print(atom.symbol, atom.position[0], atom.position[1], atom.position[2])

    # Optional (1 - spec, 0 - not spec)
    try:
        print(1, box.get_potential_energy())
    except RuntimeError:
        print(0)

    if 'config_type' in box.info:
        print(1, box.info['config_type'])
    else:
        print(0)

    varr = None
    if 'virial' in box.info:
        varr = box.info['virial']
    elif 'virials' in box.info:
        varr = box.info['virials']
    else:
        print(0)

    if varr is not None:
        print(1)
        for i in range(3):
            for j in range(3):
                print(varr[i][j])

    if 'momenta' in box.arrays:
        print(1)
        for velocity in box.get_velocities():
            print(velocity[0], velocity[1], velocity[2])
    else:
        print(0)

    farr = None
    if 'forces' in box.arrays:
        farr = box.arrays['forces']
    elif 'force' in box.arrays:
        farr = box.arrays['force']
    else:
        print(0)

    if farr is not None:
        print(1)
        for force in farr:
            print(force[0], force[1], force[2])
)";
        string placeholder = "__filename__";
        pythonScript.replace(pythonScript.find(placeholder), placeholder.size(), xyzFileName);
        return executePythonScript(pythonScript);
    }

    vector<AtomicStructure> readXyz(string_view fileName) {

        string xyzMod = xyzToJgapInput(fileName);
        istringstream inStream(xyzMod);

        size_t nBoxes;
        inStream >> nBoxes;

        vector<AtomicStructure> result(nBoxes);
        for (size_t boxIndex = 0; boxIndex < nBoxes; boxIndex++) {

            array<Vector3, 3> cell{};

            for (int i = 0; i < 3; i++) {
                    inStream >> cell[i].x >> cell[i].y >> cell[i].z;
            }

            result[boxIndex].latticeVectors = cell;

            size_t nAtoms;
            inStream >> nAtoms;

            vector<AtomData> atoms;
            for (size_t atomIndex = 0; atomIndex < nAtoms; atomIndex++) {
                string species;
                double x, y, z;
                inStream >> species >> x >> y >> z;
                atoms.push_back({
                    .position={x, y, z},
                    .species=species
                });
            }

            result[boxIndex].atoms = atoms;

            bool hasEnergy = false;
            inStream >> hasEnergy;

            if (hasEnergy) {
                double energy;
                inStream >> energy;
                result[boxIndex].energy = energy;
            }

            bool hasConfigType = false;
            inStream >> hasConfigType;

            if (hasConfigType) {
                string configType;
                inStream >> configType;
                result[boxIndex].configType = configType;
            }

            bool hasVirials = false;
            inStream >> hasVirials;
            if (hasVirials) {
                array<array<double, 3>, 3> virials{};
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        inStream >> virials[i][j];
                    }
                }
                result[boxIndex].virials = virials;
            }

            bool hasVelocities = false;
            inStream >> hasVelocities;
            if (hasVelocities) {
                for (size_t i = 0; i < nAtoms; i++) {
                    double x, y, z;
                    inStream >> x >> y >> z;
                    result[boxIndex].atoms[i].velocity = {x, y, z};
                }
            }

            bool hasForces = false;
            inStream >> hasForces;
            if (hasForces) {
                for (size_t i = 0; i < nAtoms; i++) {
                    double x, y, z;
                    inStream >> x >> y >> z;
                    result[boxIndex].atoms[i].force = {x, y, z};
                }
            }
        }

        return result;
    }

    // TODO: separate class and play around
    Vector3 toInvariantTriplet(const pair<double, double> &distanceToNodes, double distanceBetweenNodes) {
        return {
            distanceToNodes.first + distanceToNodes.second,
            pow(distanceToNodes.first - distanceToNodes.second, 2),
            distanceBetweenNodes
        };
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
}
