#ifndef JGAP_REGULARIZATIONRULES_HPP
#define JGAP_REGULARIZATIONRULES_HPP

#include "Regularization.hpp"
#include "core/atomic/neighbours/NeighbourList.hpp"
#include <memory>

namespace jgap {
    class RegularizationRules {
    public:
        virtual ~RegularizationRules() = default;
        virtual void fillSigmas(Regularization& sigmas, const Atoms& atoms) const = 0;
        virtual RegularizationRules* clone() const = 0;
    };
}

#endif