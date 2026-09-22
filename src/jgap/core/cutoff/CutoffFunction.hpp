#ifndef JGAP_CUTOFFFUNCTION_HPP
#define JGAP_CUTOFFFUNCTION_HPP

#include <string>
#include <tuple>
#include "../ValuePtr.hpp"

namespace jgap {

    class CutoffFunction {
    public:
        virtual ~CutoffFunction() = default;

        /// @note In derived classes, overriding evaluate by calling CutoffFunction::evaluate(r)
        /// (e.g. `double evaluate(double r) const override { return CutoffFunction::evaluate(r); }`)
        /// forces devirtualization of evaluateAndDifferentiate for compiler optimizations.
        virtual double evaluate(double r) const {
            return std::get<0>(evaluateAndDifferentiate(r));
        }

        virtual std::tuple<double, double> evaluateAndDifferentiate(double r) const = 0;

        virtual double getCutoff() const = 0;

        virtual CutoffFunction* clone() const = 0;
    };

    static_assert(Cloneable<CutoffFunction>);
}

#endif