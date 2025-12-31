#ifndef JGAP_TABPAIRFUNCTION_HPP
#define JGAP_TABPAIRFUNCTION_HPP

#include "EamPairFunction.hpp"
#include "../../../../data/tabulation/Grid1d.hpp"
#include "utils/BSplineTools.hpp"
#include "utils/Utils.hpp"

namespace jgap {
    class TabPairFunction : public EamPairFunction {
    public:
        TabPairFunction(const std::shared_ptr<Grid1d> &spline_grid) : spline_grid_(spline_grid) {}

        double evaluate(double distance) override {
            return BSplineTools::interpolate(*spline_grid_, distance).value;
        }
        double differentiate(double distance) override {
            return BSplineTools::interpolate(*spline_grid_, distance).gradient;
        }

    private:
        std::shared_ptr<Grid1d> spline_grid_;
    };
}

#endif