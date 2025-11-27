#ifndef JGAP_TWOBODYTABKERNEL_HPP
#define JGAP_TWOBODYTABKERNEL_HPP

#include "TwoBodyKernel.hpp"
#include "data/Grid1d.hpp"
#include "data/TabulationData.hpp"
#include "utils/BSplineTools.hpp"
#include "utils/Utils.hpp"

namespace jgap {
    class TabKernel2b : public TwoBodyKernel {
    public:
        TabKernel2b(const SpeciesPair &idSpecies, const std::shared_ptr<Grid1d> &splineGrid)
            : _idSpecies(idSpecies), _splineGrid(splineGrid) {
            coefficient = 1.0;
        }

        std::string getType() override { IO_NOT_INTENDED(TabKernel2b.getType); }
        nlohmann::json serialize() override { IO_NOT_INTENDED(TabKernel2b.serialize); }
        double crossCovariance(const std::shared_ptr<IKernel> &kernel) override {
            throw std::logic_error("TabKernel2b not intended for fitting | crossCovariance()");
        }
        SpeciesPair getFilter() override { return _idSpecies; };

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
