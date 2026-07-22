#ifndef JGAP_TWOBODYTGCOMPONENT_HPP
#define JGAP_TWOBODYTGCOMPONENT_HPP

#include "TabGapComponent.hpp"
#include "core/ValuePtr.hpp"
#include "core/atomic/iteration/Cluster2Expansion.hpp"
#include "core/atomic/species/composition/Species2Sorted.hpp"
#include "core/splines/Spline.hpp"

namespace jgap {
    class TwoBodyTGComponent : public TabGapComponent {
    public:
        TwoBodyTGComponent(Species2Sorted species_pair, ValuePtr<Spline<1>> spline);

        AtomicQuantity energy(const NeighbourLists& nl) const override;

        Cutoffs getCutoffs() const override;

        void tabulate(TabulationData& tables) const override;

        std::set<Species> getAllSpecies() const override { return {species_pair.nodes[0], species_pair.nodes[1]}; }

        TwoBodyTGComponent* clone() const override { return new TwoBodyTGComponent(*this); }

        const auto& getSpeciesPair() const { return species_pair; }

        const auto& getSpline() const { return spline; }

    private:
        Species2Sorted species_pair;
        ValuePtr<Spline<1>> spline;
        Cluster2Expansion expansion;
    };
}

#endif
