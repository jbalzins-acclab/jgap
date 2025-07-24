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

using namespace std;

namespace jgap {

    class Descriptor/*(Manager)*/ {
    public:

        virtual ~Descriptor() = default;

        // Sparsification strategy to constructor
        virtual void setSparsePoints(const vector<AtomicStructure> &fromData) = 0;
        virtual size_t nSparsePoints() = 0;

        virtual vector<Covariance> covariate(const AtomicStructure &atomicStructure) = 0;
        virtual vector<pair<size_t/*sparse point id*/, shared_ptr<MatrixBlock>>> selfCovariate() = 0;

        virtual nlohmann::json serialize() = 0;
        virtual string getType() = 0;
        virtual double getCutoff() = 0;

        virtual TabulationData tabulate(const TabulationParams &params) = 0;

        void setCoefficients(const vector<double>& c);

        PotentialPrediction predict(const AtomicStructure &atomicStructure);

    protected:
        vector<double> _coefficients = {};
    };
}

#endif
