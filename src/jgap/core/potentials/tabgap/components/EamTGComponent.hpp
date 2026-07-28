#ifndef JGAP_EAMTGCOMPONENT_HPP
#define JGAP_EAMTGCOMPONENT_HPP
#include "TabGapComponent.hpp"
#include "jgap/core/splines/HermiteCubicSpline.hpp"
#include "jgap/core/transform/manybody/TwoBodySum.hpp"

namespace jgap {
    class EamTGComponent : public TabGapComponent {
    public:
        EamTGComponent(ManyBodyGrids2<1, 1> grids_for_element);

        AtomicQuantity energy(const NeighbourLists& nl) const override;

        Cutoffs getCutoffs() const override;

        void tabulate(TabulationData& tables) const override;

        std::set<Species> getAllSpecies() const override;

        EamTGComponent* clone() const override { return new EamTGComponent(*this); }

        Species getCentralSpecies() const { return spline_density_aggregator.getCentralSpecies(); }

        auto getSplineNBodyAggregator() const { return spline_density_aggregator; }

        auto getEnergySpline() const { return energy_spline; }

    private:
        TwoBodySum<1> spline_density_aggregator;
        HermiteCubicSpline energy_spline;
    };
}

#endif
