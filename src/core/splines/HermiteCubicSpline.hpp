#ifndef JGAP_HERMITECUBICSPLINE_HPP
#define JGAP_HERMITECUBICSPLINE_HPP

#include <vector>
#include "core/Real.hpp"
#include "Spline.hpp"
#include "Table.hpp"

namespace jgap {
    class HermiteCubicSpline : public Spline<1> {
    public:
        explicit HermiteCubicSpline(const Table<1>& table);

        InterpolationResults<1> interpolate(std::array<Real, 1> pos) const override;
        std::array<Real, 1> getCutoff() const override;

    private:
        void init(const Table<1>& table);
        size_t findInterval(Real r) const;

        Table<1> table;
        std::vector<Real> b;
        std::vector<Real> c;
        std::vector<Real> d;
    };
}

#endif