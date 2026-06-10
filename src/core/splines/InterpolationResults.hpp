#ifndef JGAP_INTERPOLATIONRESULTS_HPP
#define JGAP_INTERPOLATIONRESULTS_HPP
#include <array>

#include "core/Real.hpp"

namespace jgap {
    template<size_t Dim>
    struct InterpolationResults {
        Real value;
        std::array<Real, Dim> gradient;
    };
}

#endif
