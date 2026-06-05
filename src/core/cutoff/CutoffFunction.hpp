#ifndef JGAP_CUTOFFFUNCTION_HPP
#define JGAP_CUTOFFFUNCTION_HPP

#include <string>
#include <tuple>
#include "core/Real.hpp"

namespace jgap {

    class CutoffFunction {
    public:
        virtual ~CutoffFunction() = default;

        virtual Real evaluate(Real r) const = 0;
        virtual std::tuple<Real, Real> evaluateAndDifferentiate(Real r) const = 0;

        virtual Real getCutoff() const = 0;
    };
}

#endif