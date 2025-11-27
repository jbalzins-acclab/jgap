#ifndef JGAP_FIT_HPP
#define JGAP_FIT_HPP

#include "data/Vector3.hpp"
#include "core/potentials/Potential.hpp"

#include <memory>
#include <core/tabulation/SimpleTabulation.hpp>


namespace jgap {
    class Fit {
    public:
        virtual ~Fit() = default;
        virtual shared_ptr<Potential> fit(const vector<AtomicStructure>& trainingData) = 0;
    };
}

#endif
