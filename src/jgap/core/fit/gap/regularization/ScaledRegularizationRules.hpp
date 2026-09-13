#ifndef JGAP_SCALEDREGULARIZATIONRULES_HPP
#define JGAP_SCALEDREGULARIZATIONRULES_HPP

#include <memory>
#include "PerConfigTypeSigmas.hpp"
#include "RegularizationRules.hpp"

namespace jgap {

    /// \brief Regularization rules that dynamically scale target sigmas (uncertainties)
    /// with the magnitude of the maximum atomic force in each structure.
    ///
    /// Given maximum force magnitude F_max = max_i ||F_i||:
    ///   scale = max(min_scale, 1.0 + force_scale * F_max)
    ///   sigma_scaled = sigma_base * scale
    class ScaledRegularizationRules : public RegularizationRules {
    public:
        explicit ScaledRegularizationRules(
            std::shared_ptr<RegularizationRules> base_rules,
            double force_scale = 1.0,
            double min_scale = 1.0
        );

        explicit ScaledRegularizationRules(
            const RegularizationRules& base_rules,
            double force_scale = 1.0,
            double min_scale = 1.0
        );

        explicit ScaledRegularizationRules(
            PerConfigTypeSigmas base_sigmas,
            double force_scale = 1.0,
            double min_scale = 1.0
        );

        explicit ScaledRegularizationRules(
            double energy_sigma_per_atom = 0.001,
            double force_component_sigma = 0.05,
            double virials_iso_sigma_per_atom = 0.1,
            double virials_aniso_sigmas_per_atom = 0.02,
            double force_scale = 1.0,
            double min_scale = 1.0
        );

        Regularization determine(const Atoms& atoms) const override;

        ScaledRegularizationRules* clone() const override {
            return new ScaledRegularizationRules(*this);
        }

        double getForceScale() const { return force_scale; }
        double getMinScale() const { return min_scale; }
        const RegularizationRules* getBaseRules() const { return base_rules.get(); }

    private:
        std::shared_ptr<RegularizationRules> base_rules;
        double force_scale;
        double min_scale;
    };

}

#endif
