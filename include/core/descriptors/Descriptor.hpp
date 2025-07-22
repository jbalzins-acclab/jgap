#ifndef DESCRIPTOR_HPP
#define DESCRIPTOR_HPP

#include "data/BasicDataTypes.hpp"
#include "io/log/CurrentLogger.hpp"
#include "memory/MatrixBlock.hpp"

#include <memory>
#include <utility>
#include <vector>
#include <nlohmann/json_fwd.hpp>

#include "data/TabulationData.hpp"
#include "utils/Utils.hpp"

using namespace std;

namespace jgap {

    class Descriptor {
    public:

        virtual ~Descriptor() = default;

        // Sparsification strategy to constructor
        virtual void setSparsePoints(const vector<AtomicStructure> &fromData) = 0;
        virtual double getCutoff() = 0;
        virtual size_t nSparsePoints() = 0;

        virtual vector<Covariance> covariate(const AtomicStructure &atomicStructure) = 0;
        virtual vector<pair<size_t/*sparse point id*/, shared_ptr<MatrixBlock>>> selfCovariate() = 0;

        virtual nlohmann::json serialize() = 0;
        virtual string getType() = 0;

        virtual TabulationData tabulate(const TabulationParams &params) = 0;

        void setCoefficients(const vector<double>& c) {
            if (c.size() != nSparsePoints()) {
                CurrentLogger::get()->error("Number of coefficients doesn't match the number of sparse points");
            }
            _coefficients = c;
        }

        PotentialPrediction predict(const AtomicStructure &atomicStructure) {
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

    protected:
        vector<double> _coefficients = {};
    };
}

#endif
