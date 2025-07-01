//
// Created by Jegors Balzins on 23.6.2025.
//

#ifndef SIGMARULES_HPP
#define SIGMARULES_HPP
#include "data/BasicDataTypes.hpp"

namespace jgap {
    class SigmaRules {
        public:
        virtual ~SigmaRules() = default;
        virtual void fillSigmas(AtomicStructure& structure) = 0;
    };
}

#endif //SIGMARULES_HPP
