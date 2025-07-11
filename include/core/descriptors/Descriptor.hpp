#ifndef DESCRIPTOR_HPP
#define DESCRIPTOR_HPP

#include "data/BasicDataTypes.hpp"
#include "io/log/CurrentLogger.hpp"
#include "memory/MatrixBlock.hpp"

#include <memory>
#include <utility>
#include <vector>
#include <nlohmann/json_fwd.hpp>

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
            vector forces(atomicStructure.atoms.size(), Vector3(0.0, 0.0, 0.0));

            for (size_t i = 0; i < _coefficients.size(); i++) {
                energy += _coefficients[i] * covariance[i].total;
                for (size_t j = 0; j < atomicStructure.atoms.size(); j++) {
                    forces[j] = forces[j] - covariance[i].derivatives[j] * _coefficients[i];
                }
            }

            return {
                energy,
                forces,
                {}
            };
        }

    protected:
        vector<double> _coefficients = {};
    };
}

#endif
