#ifndef JGAP_TABKERNELEAM_HPP
#define JGAP_TABKERNELEAM_HPP

#include "EamKernel.hpp"
#include "data/Grid1d.hpp"
#include "utils/BSplineTools.hpp"
#include "utils/Utils.hpp"

namespace jgap {
    class TabKernelEam : public EamKernel {
    public:
        TabKernelEam(const Species &idSpecies, const shared_ptr<Grid1d> &grid)
            : _idSpecies(idSpecies), _splineGrid(grid) {
            coefficient = 1.0;
        }

        string getType() override { IO_NOT_INTENDED(TabKernelEam.getType); }
        nlohmann::json serialize() override { IO_NOT_INTENDED(TabKernelEam.serialize); }
        double crossCovariance(const shared_ptr<IKernel> &kernel) override {
            throw logic_error("TabKernel3b not intended for fitting | crossCovariance()");
        }

        string getFilter() override { return _idSpecies; }
        pair<double, double> getDensityRange() override {
            return {_splineGrid->origin, _splineGrid->cutoff()};
        }

    protected:
        double valueInternal(const double &density) const override {
            return BSplineTools::interpolate(*_splineGrid, density).value;
        }
        double derivativeInternal(const double &density) const override {
            return BSplineTools::interpolate(*_splineGrid, density).gradient;
        }

    private:
        Species _idSpecies;
        shared_ptr<Grid1d> _splineGrid;
    };
}

#endif
