#ifndef JGAP_CUTOFFFUNCTION_HPP
#define JGAP_CUTOFFFUNCTION_HPP

#include <string>
#include <tuple>
#include "../ValuePtr.hpp"
#include "jgap/core/Real.hpp"

namespace jgap {

    class CutoffFunction {
    public:
        virtual ~CutoffFunction() = default;

        /// @note In derived classes, overriding evaluate by calling CutoffFunction::evaluate(r)
        /// (e.g. `Real evaluate(Real r) const override { return CutoffFunction::evaluate(r); }`)
        /// forces devirtualization of evaluateAndDifferentiate for compiler optimizations.
        virtual Real evaluate(Real r) const {
            return std::get<0>(evaluateAndDifferentiate(r));
        }

        virtual std::tuple<Real, Real> evaluateAndDifferentiate(Real r) const = 0;

        virtual Real getCutoff() const = 0;

        virtual CutoffFunction* clone() const = 0;
    };

    static_assert(Cloneable<CutoffFunction>);
}

#endif