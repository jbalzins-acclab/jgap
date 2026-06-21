#ifndef JGAP_EAMTGCOMPONENT_HPP
#define JGAP_EAMTGCOMPONENT_HPP
#include "TabGapComponent.hpp"
#include "core/splines/HermiteCubicSpline.hpp"
#include "core/transform/3b/CosTransformation.hpp"
#include "core/transform/aggregated/TransformationAggregator.hpp"
#include "core/transform/aggregated/TransformationAggregatorImpl.hpp"

namespace jgap {
    class EamTGComponent : public TabGapComponent {
    public:
        EamTGComponent(ManyBodyGrids<2> grids_for_element);

        AtomicQuantity energy(const NeighbourList &nl) const override;

        Cutoffs getCutoffs() const override;

        void tabulate(TabulationData &tables) const override;

        std::set<Species> getAllSpecies() const override;

        EamTGComponent* clone() const override {
            return new EamTGComponent(*this);
        }

        auto getSplineTransformationAggregator() const {
            return spline_density_aggregator;
        }

        auto getEnergySpline() const {
            return energy_spline;
        }

    private:
        TransformationAggregatorImpl<1, 2> spline_density_aggregator;
        HermiteCubicSpline energy_spline;
    };
}

#endif
