#ifndef JGAP_CUBICBSPLINE_HPP
#define JGAP_CUBICBSPLINE_HPP

#include <vector>

#include "InterpolationResults.hpp"
#include "Spline.hpp"
#include "Grid.hpp"
#include "../ValuePtr.hpp"

namespace jgap {
    class CubicBSpline : public Spline<1> {
    public:
        static std::vector<double> toSplineCoefficients(const std::vector<double>& original, double spacing);

        static ValuePtr<Spline> fit(const Grid<1>& values);

        explicit CubicBSpline(const Grid<1>& coefficients);

        InterpolationResults<1> interpolate(std::array<double, 1> pos) const override;
        std::array<double, 1> getCutoff() const override;

        CubicBSpline* clone() const override {
            return new CubicBSpline(*this);
        }

        const Grid<1>& getCoefficients() const {
            return coefficients;
        }

    private:
        Grid<1> coefficients;
    };
}

#endif