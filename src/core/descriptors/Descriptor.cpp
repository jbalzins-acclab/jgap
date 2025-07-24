#include "core/descriptors/Descriptor.hpp"
#include "data/BasicDataTypes.hpp"

using namespace std;

namespace jgap {

    void Descriptor::setCoefficients(const vector<double> &c) {
        if (c.size() != nSparsePoints()) {
            CurrentLogger::get()->error("Number of coefficients doesn't match the number of sparse points");
        }
        _coefficients = c;
    }

    PotentialPrediction Descriptor::predict(const AtomicStructure &atomicStructure) {
        if (_coefficients.size() != nSparsePoints()) {
            CurrentLogger::get()->error("Coefficients not set");
        }

        const auto covariance = covariate(atomicStructure);

        double energy = 0.0;
        vector forces(atomicStructure.size(), Vector3(0.0, 0.0, 0.0));
        array<Vector3, 3> virials{};

        for (size_t i = 0; i < _coefficients.size(); i++) {
            energy += _coefficients[i] * covariance[i].total;

            for (size_t j = 0; j < atomicStructure.size(); j++) {
                forces[j] += covariance[i].forces[j] * _coefficients[i];
            }

            virials[0] += covariance[i].virials[0] * _coefficients[i];
            virials[1] += covariance[i].virials[1] * _coefficients[i];
            virials[2] += covariance[i].virials[2] * _coefficients[i];
        }

        return {
            energy,
            forces,
            virials
        };
    }
}