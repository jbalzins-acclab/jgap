#ifndef JGAP_SPLINE_HPP
#define JGAP_SPLINE_HPP

#include <memory>

#include "InterpolationResults.hpp"

namespace jgap {

    template<size_t Dim>
    class Spline {
    public:
        virtual ~Spline() = default;
        virtual InterpolationResults<Dim> interpolate(std::array<Real, Dim> pos) const = 0;
        virtual std::array<Real, Dim> getCutoff() const = 0;
        virtual std::unique_ptr<Spline<Dim>> clone() const = 0;
    };
}

#endif
