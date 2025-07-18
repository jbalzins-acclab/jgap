#ifndef JGAPFIT_HPP
#define JGAPFIT_HPP

#include "data/BasicDataTypes.hpp"
#include "io/log/CurrentLogger.hpp"
#include "core/potentials/Potential.hpp"
#include "core/descriptors/Descriptor.hpp"
#include "core/potentials/JgapPotential.hpp"
#include "core/matrices/sigmas/SigmaRules.hpp"

#include <Eigen/Dense>

#include <memory>

using namespace std;

namespace jgap {
    class Fit {
    public:
        virtual ~Fit() = default;
        virtual shared_ptr<Potential> fit(const vector<AtomicStructure>& trainingData) = 0;
        virtual string getType() = 0; // for logging purposes
    };
}
#endif //JGAPFIT_HPP
