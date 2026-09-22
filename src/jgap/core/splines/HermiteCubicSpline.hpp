#ifndef JGAP_HERMITECUBICSPLINE_HPP
#define JGAP_HERMITECUBICSPLINE_HPP

#include <vector>
#include "Grid.hpp"
#include "Spline.hpp"

namespace jgap {
    class HermiteCubicSpline : public Spline<1> {
    public:
        explicit HermiteCubicSpline(const Grid<1>& table);

        InterpolationResults<1> interpolate(std::array<double, 1> pos) const override;
        std::array<double, 1> getCutoff() const override;

        HermiteCubicSpline* clone() const override { return new HermiteCubicSpline(*this); }

        const Grid<1>& getTable() const { return table; }

    private:
        Grid<1> table;
        std::vector<double> b;
        std::vector<double> c;
        std::vector<double> d;

        void init(const Grid<1>& table);
        size_t findInterval(double r) const;
    };
}

#endif
