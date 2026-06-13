#ifndef JGAP_CUBICBSPLINE3D_HPP
#define JGAP_CUBICBSPLINE3D_HPP

#include "InterpolationResults.hpp"
#include "Spline.hpp"
#include "Grid.hpp"
#include "core/Real.hpp"
#include "core/atomic/geometry/Vector3.hpp"

namespace jgap {
    class CubicBSpline3D : public Spline<3> {
    public:
        explicit CubicBSpline3D(const Grid<3>& coefficients);

        static CubicBSpline3D fit(const Grid<3>& values);

        InterpolationResults<3> interpolate(std::array<Real, 3> pos) const override;
        std::array<Real, 3> getCutoff() const override;

        std::unique_ptr<Spline<3>> clone() const override {
            return std::make_unique<CubicBSpline3D>(*this);
        }

        const Grid<3>& getCoefficients() const {
            return coefficients;
        }

    private:
        Grid<3> coefficients;
    };
}

#endif