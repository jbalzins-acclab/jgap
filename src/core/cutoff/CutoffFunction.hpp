#ifndef JGAP_CUTOFFFUNCTION_HPP
#define JGAP_CUTOFFFUNCTION_HPP

#include <string>
#include <tuple>
#include "core/Real.hpp"
#include "../ValuePtr.hpp"

namespace jgap {

    class CutoffFunction {
    public:
        // Out-of-line key function: forces a single externally-linked vtable/type_info so dynamic_cast
        // works across the jgap_lib boundary (libc++abi compares type_info by address on Apple).
        // Defined, together with the other serializable types' anchors, in SerializationTypeAnchors.cpp.
        virtual ~CutoffFunction();

        virtual Real evaluate(Real r) const = 0;
        virtual std::tuple<Real, Real> evaluateAndDifferentiate(Real r) const = 0;

        virtual Real getCutoff() const = 0;

        virtual std::unique_ptr<CutoffFunction> clone() const = 0;
    };

    static_assert(Cloneable<CutoffFunction>);
}

#endif