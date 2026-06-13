#ifndef JGAP_SPLINEPAIRPOTENTIAL_HPP
#define JGAP_SPLINEPAIRPOTENTIAL_HPP

#include "../../splines/NaturalCubicSpline.hpp"
#include "core/atomic/species/SpeciesSet.hpp"
#include "core/potentials/Potential.hpp"
#include "io/Serializable.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class SplinePairPotential : public Potential {
    public:
        SplinePairPotential() = default;
        ~SplinePairPotential() override = default;

        void extend(Species species1, Species species2, const std::vector<Real>& r, const std::vector<Real>& energies);

        AtomicQuantity calculateEnergy(const Atoms &atoms) const override;

        Cutoffs getCutoffs() const override;

        void tabulate(TabulationData &tables) const override;

        std::unique_ptr<Potential> clone() const override {
            return std::make_unique<SplinePairPotential>(*this);
        }

    private:
        std::map<SpeciesSet<2, Symmetric>, NaturalCubicSpline> per_species_interpolators;
    };
}

#endif
