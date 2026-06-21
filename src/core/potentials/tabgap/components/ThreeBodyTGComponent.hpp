#ifndef JGAP_THREEBODYTGCOMPONENT_HPP
#define JGAP_THREEBODYTGCOMPONENT_HPP
#include "TabGapComponent.hpp"
#include "core/splines/CubicBSpline3D.hpp"
#include "core/transform/3b/CosTransformation.hpp"

namespace jgap {

    class ThreeBodyTGComponent : public TabGapComponent {
    public:

        ThreeBodyTGComponent(const SpeciesSet<3, HasCentralAtom> &species, const CubicBSpline3D &spline);

        AtomicQuantity energy(const NeighbourList &nl) const override;

        Cutoffs getCutoffs() const override;

        void tabulate(TabulationData &tables) const override;

        std::set<Species> getAllSpecies() const override;

        ThreeBodyTGComponent* clone() const override {
            return new ThreeBodyTGComponent(*this);
        }

        auto getSpeciesTriplet() const {
            return species;
        }

        const auto& getSpline() const {
            return spline;
        }

    private:
        SpeciesSet<3, HasCentralAtom> species;
        CubicBSpline3D spline;

        CosTransformation transformation{};
    };
}

#endif
