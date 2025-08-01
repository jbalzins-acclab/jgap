#ifndef REGULARIZATIONRULES_HPP
#define REGULARIZATIONRULES_HPP

#include "data/BasicDataTypes.hpp"

namespace jgap {
    class RegularizationRules {
    public:
        virtual ~RegularizationRules() = default;
        virtual void fillSigmas(AtomicStructure& structure) = 0;
    };
}

#endif
