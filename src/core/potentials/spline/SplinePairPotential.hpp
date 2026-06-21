#ifndef JGAP_SPLINEPAIRPOTENTIAL_HPP
#define JGAP_SPLINEPAIRPOTENTIAL_HPP

#include "../../splines/NaturalCubicSpline.hpp"
#include "core/atomic/species/SpeciesSet.hpp"
#include "core/potentials/Potential.hpp"
#include "../../../serialization/SerializationRegistry.hpp"

namespace jgap {
    class SplinePairPotential : public Potential {
    public:
        SplinePairPotential() = default;
        ~SplinePairPotential() override = default;

        void extend(Species species1, Species species2, const std::vector<Real>& r, const std::vector<Real>& energies);

        AtomicQuantity calculateEnergy(const Atoms &atoms) const override;

        Cutoffs getCutoffs() const override;

        void fillTables(TabulationData &tables) const override;

        SplinePairPotential* clone() const override {
            return new SplinePairPotential(*this);
        }

        const auto& getInterpolators() const {
            return per_species_interpolators;
        }

    private:
        std::map<SpeciesSet<2, FullSymmetry>, NaturalCubicSpline> per_species_interpolators;
    };
}

#endif