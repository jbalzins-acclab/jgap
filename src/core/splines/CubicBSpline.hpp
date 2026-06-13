#ifndef JGAP_CUBICBSPLINE_HPP
#define JGAP_CUBICBSPLINE_HPP

#include <vector>

#include "InterpolationResults.hpp"
#include "Spline.hpp"
#include "Grid.hpp"
#include "core/Real.hpp"
#include "utils/ValuePtr.hpp"

namespace jgap {
    class CubicBSpline : public Spline<1> {
    public:
        static std::vector<Real> toSplineCoefficients(const std::vector<Real>& original, Real spacing);

        static ValuePtr<Spline> fit(const Grid<1>& values);

        explicit CubicBSpline(const Grid<1>& coefficients);

        InterpolationResults<1> interpolate(std::array<Real, 1> pos) const override;
        std::array<Real, 1> getCutoff() const override;

        std::unique_ptr<Spline<1>> clone() const override {
            return std::make_unique<CubicBSpline>(*this);
        }

        const Grid<1>& getCoefficients() const {
            return coefficients;
        }

    private:
        Grid<1> coefficients;
    };
}

#endif