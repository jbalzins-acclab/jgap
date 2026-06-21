#ifndef JGAP_TWOBODYTGCOMPONENT_HPP
#define JGAP_TWOBODYTGCOMPONENT_HPP

#include "TabGapComponent.hpp"
#include "core/splines/CubicBSpline.hpp"
#include "core/splines/Spline.hpp"

namespace jgap {
    class TwoBodyTGComponent : public TabGapComponent {
    public:
        TwoBodyTGComponent(SpeciesSet<2, FullSymmetry> species_pair, ValuePtr<Spline<1>> spline);

        AtomicQuantity energy(const NeighbourList &nl) const override;

        Cutoffs getCutoffs() const override;

        void tabulate(TabulationData &tables) const override;

        std::set<Species> getAllSpecies() const override {
            return {species_pair.getNodes()[0], species_pair.getNodes()[1]};
        }

        TwoBodyTGComponent* clone() const override {
            return new TwoBodyTGComponent(*this);
        }

        const auto& getSpeciesPair() const {
            return species_pair;
        }

        const auto& getSpline() const {
            return spline;
        }

    private:
        SpeciesSet<2, FullSymmetry> species_pair;
        ValuePtr<Spline<1>> spline;
    };
}

#endif
