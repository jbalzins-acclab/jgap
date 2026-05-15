#ifndef JGAP_SIMPLEREGULARIZATIONRULES_HPP
#define JGAP_SIMPLEREGULARIZATIONRULES_HPP

#include "Regularization.hpp"
#include "RegularizationRules.hpp"
#include "core/atomic/neighbours/NeighbourList.hpp"

namespace jgap {
    class SimpleRegularizationRules : public RegularizationRules {
    public:
        SimpleRegularizationRules(Real energy_sigma_per_atom = 0.001,
                                  Real force_component_sigma = 0.05,
                                  Real virials_iso_sigma_per_atom = 0.1,
                                  Real virials_aniso_sigmas_per_atom = 0.02,
                                  Real liquid_multiplier = 5.0,
                                  Real short_range_multiplier = 5.0);

        void fillSigmas(Regularization& sigmas, const Atoms& atoms, const NeighbourList& neighbour_list) override;

    private:
        Real default_energy_per_atom;
        Vector3 default_force;
        Virials default_virials_per_atom;

        Real liquid_multiplier;
        Real short_range_multiplier;
    };
}
#endif
