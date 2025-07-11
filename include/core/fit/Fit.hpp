
//
// Created by Jegors Balzins on 9.7.2025.
//

#ifndef JGAPFIT_HPP
#define JGAPFIT_HPP

#include "data/BasicDataTypes.hpp"
#include "io/log/CurrentLogger.hpp"
#include "../potentials/Potential.hpp"
#include "core/descriptors/Descriptor.hpp"
#include "../potentials/JgapPotential.hpp"
#include "../matrices/sigmas/SigmaRules.hpp"

#include <Eigen/Dense>

#include <memory>
#include <concepts>

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
