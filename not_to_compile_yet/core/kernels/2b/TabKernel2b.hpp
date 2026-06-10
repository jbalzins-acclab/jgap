#ifndef JGAP_TWOBODYTABKERNEL_HPP
#define JGAP_TWOBODYTABKERNEL_HPP

#include <utility>

#include "../descriptors/2b/TwoBodyKernel.hpp"
#include "../../../../src/core/splines/Grid1d.hpp"
#include "../../../../src/core/tabulation/TabulationData.hpp"
#include "utils/BSplineTools.hpp"
#include "utils/Utils.hpp"

namespace jgap {
    class TabKernel2b : public TwoBodyKernel {
    public:
        TabKernel2b(SpeciesPair idSpecies, const std::shared_ptr<Grid1d> &splineGrid)
            : _idSpecies(std::move(idSpecies)), _splineGrid(splineGrid) {
            coefficient = 1.0;
        }

        double crossCovariance(const std::shared_ptr<IKernel> &kernel) override {
            throw std::logic_error("TabKernel2b not intended for fitting | crossCovariance()");
        }
        SpeciesPair getFilter() override { return _idSpecies; }

    protected:
        double valueInternal(const double &r) const override {
            return BSplineTools::interpolate(*_splineGrid, r).value;
        }
        double derivativeInternal(const double &changingR) const override {
            return BSplineTools::interpolate(*_splineGrid, changingR).gradient;
        }

    private:
        SpeciesPair _idSpecies;
        std::shared_ptr<Grid1d> _splineGrid;
    };
}

#endif
