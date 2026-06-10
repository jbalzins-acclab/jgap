#ifndef JGAP_CUBICBSPLINE_HPP
#define JGAP_CUBICBSPLINE_HPP

#include <vector>

#include "InterpolationResults.hpp"
#include "Spline.hpp"
#include "Table.hpp"
#include "core/Real.hpp"

namespace jgap {
    class CubicBSpline : public Spline<1> {
    public:
        static std::vector<Real> toSplineCoefficients(const std::vector<Real>& original, Real spacing);

        static CubicBSpline fit(const Table<1>& values);

        explicit CubicBSpline(const Table<1>& coefficients);

        InterpolationResults<1> interpolate(std::array<Real, 1> pos) const override;
        std::array<Real, 1> getCutoff() const override;

    private:
        Table<1> coefficients;
    };
}

#endif