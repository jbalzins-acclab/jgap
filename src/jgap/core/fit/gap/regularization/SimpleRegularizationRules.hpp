#ifndef JGAP_SIMPLEREGULARIZATIONRULES_HPP
#define JGAP_SIMPLEREGULARIZATIONRULES_HPP

#include "PerConfigTypeSigmas.hpp"
#include "RegularizationRules.hpp"

namespace jgap {
    class SimpleRegularizationRules : public RegularizationRules {
    public:
        SimpleRegularizationRules(
            double energy_sigma_per_atom = 0.001,
            double force_component_sigma = 0.05,
            double virials_iso_sigma_per_atom = 0.1,
            double virials_aniso_sigmas_per_atom = 0.02,
            double liquid_multiplier = 5.0,
            double short_range_multiplier = 5.0
        );

        Regularization determine(const Atoms& atoms) const override;

        SimpleRegularizationRules* clone() const override {
            return new SimpleRegularizationRules(*this);
        }

    private:
        PerConfigTypeSigmas defaults;

        double liquid_multiplier;
        double short_range_multiplier;
    };
}
#endif