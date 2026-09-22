#ifndef JGAP_INTERPOLATIONRESULTS_HPP
#define JGAP_INTERPOLATIONRESULTS_HPP
#include <array>


namespace jgap {
    template<size_t Dim>
    struct InterpolationResults {
        double value;
        std::array<double, Dim> gradient;
    };
}

#endif
