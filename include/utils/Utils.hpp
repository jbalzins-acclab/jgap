#ifndef UTILS_HPP
#define UTILS_HPP

#include "data/BasicDataTypes.hpp"

#include <string>
#include <string_view>
#include <vector>
#include <memory>
#include <Eigen/Dense>

using namespace std;

namespace jgap {
    string runBash(string_view command);
    string executePythonScript(string_view pythonCode);
    string xyzToJgapInput(string_view xyzFileName);
    vector<AtomicStructure> readXyz(string_view fileName);
    Vector3 toInvariantTriplet(const pair<double, double> &distanceToNodes, double distanceBetweenNodes);
    vector<string> split(const string& s, char delimiter);
    void saveArray(const vector<double>& data, const string& filename);
    vector<double> loadArray(const string& filename);
    string matrixToString(const Eigen::MatrixXd& mat);
}

#endif
