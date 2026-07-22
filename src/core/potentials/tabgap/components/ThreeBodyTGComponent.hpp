#ifndef JGAP_THREEBODYTGCOMPONENT_HPP
#define JGAP_THREEBODYTGCOMPONENT_HPP
#include "TabGapComponent.hpp"
#include "core/atomic/iteration/AtomicCluster3Expansion.hpp"
#include "core/atomic/species/composition/Species3AtomicSorted.hpp"
#include "core/splines/CubicBSpline3D.hpp"
#include "core/transform/nbody/3b/CosTransformation.hpp"

namespace jgap {

    class ThreeBodyTGComponent : public TabGapComponent {
    public:
        ThreeBodyTGComponent(const Species3AtomicSorted& species, const CubicBSpline3D& spline);

        AtomicQuantity energy(const NeighbourLists& nl) const override;

        Cutoffs getCutoffs() const override;

        void tabulate(TabulationData& tables) const override;

        std::set<Species> getAllSpecies() const override;

        ThreeBodyTGComponent* clone() const override { return new ThreeBodyTGComponent(*this); }

        auto getSpeciesTriplet() const { return species; }

        const auto& getSpline() const { return spline; }

    private:
        Species3AtomicSorted species;
        CubicBSpline3D spline;

        CosTransformation transformation{};
        AtomicCluster3Expansion expansion;
    };
}

#endif
