#ifndef JGAP_SIMPLEREGULARIZATIONRULES_HPP
#define JGAP_SIMPLEREGULARIZATIONRULES_HPP

#include "RegularizationRules.hpp"
#include "ConfigSigmas.hpp"
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

        void fillSigmas(Regularization& sigmas, const Atoms& atoms) const override;

        std::unique_ptr<RegularizationRules> clone() const override {
            return std::make_unique<SimpleRegularizationRules>(*this);
        }

    private:
        ConfigSigmas defaults;

        Real liquid_multiplier;
        Real short_range_multiplier;
    };
}
#endif