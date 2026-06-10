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

        AtomicQuantity calculateEnergy(const Atoms &atoms) override;

        Cutoffs getCutoffs() override;

    private:
        std::map<SpeciesSet<2, Symmetric>, NaturalCubicSpline> per_species_interpolators;
    };
}

#endif
