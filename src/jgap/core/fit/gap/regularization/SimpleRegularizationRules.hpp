#ifndef JGAP_SIMPLEREGULARIZATIONRULES_HPP
#define JGAP_SIMPLEREGULARIZATIONRULES_HPP

#include "PerConfigTypeSigmas.hpp"
#include "RegularizationRules.hpp"

namespace jgap {
    class SimpleRegularizationRules : public RegularizationRules {
    public:
        SimpleRegularizationRules(
            Real energy_sigma_per_atom = 0.001,
            Real force_component_sigma = 0.05,
            Real virials_iso_sigma_per_atom = 0.1,
            Real virials_aniso_sigmas_per_atom = 0.02,
            Real liquid_multiplier = 5.0,
            Real short_range_multiplier = 5.0
        );

        Regularization determine(const Atoms& atoms) const override;

        SimpleRegularizationRules* clone() const override {
            return new SimpleRegularizationRules(*this);
        }

    private:
        PerConfigTypeSigmas defaults;

        Real liquid_multiplier;
        Real short_range_multiplier;
    };
}
#endif