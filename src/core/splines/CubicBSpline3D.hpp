#ifndef JGAP_CUBICBSPLINE3D_HPP
#define JGAP_CUBICBSPLINE3D_HPP

#include "InterpolationResults.hpp"
#include "Spline.hpp"
#include "Table.hpp"
#include "core/Real.hpp"
#include "core/atomic/geometry/Vector3.hpp"

namespace jgap {
    class CubicBSpline3D : public Spline<3> {
    public:
        explicit CubicBSpline3D(const Table<3>& coefficients);

        static CubicBSpline3D fit(const Table<3>& values);

        InterpolationResults<3> interpolate(std::array<Real, 3> pos) const override;
        std::array<Real, 3> getCutoff() const override;

    private:
        Table<3> coefficients;
    };
}

#endif