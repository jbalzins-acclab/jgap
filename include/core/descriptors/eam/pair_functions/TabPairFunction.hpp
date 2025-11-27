#ifndef JGAP_TABPAIRFUNCTION_HPP
#define JGAP_TABPAIRFUNCTION_HPP

#include "EamPairFunction.hpp"
#include "data/Grid1d.hpp"
#include "utils/BSplineTools.hpp"
#include "utils/Utils.hpp"

namespace jgap {
    class TabPairFunction : public EamPairFunction {
    public:
        static constexpr string TYPE = "tab";

        TabPairFunction(shared_ptr<Grid1d> splineGrid) : _splineGrid(splineGrid) {};

        string getType() override { IO_NOT_INTENDED(TabPairFunction.getType); }
        nlohmann::json serialize() override { IO_NOT_INTENDED(TabPairFunction.serialize); }

        double evaluate(double distance) override {
            return BSplineTools::interpolate(*_splineGrid, distance).value;
        }
        double differentiate(double distance) override {
            return BSplineTools::interpolate(*_splineGrid, distance).gradient;
        }

    private:
        shared_ptr<Grid1d> _splineGrid;
    };
}

#endif