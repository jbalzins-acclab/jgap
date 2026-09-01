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
            Real force_scale = 1.0_r,
            Real min_scale = 1.0_r
        );

        explicit ScaledRegularizationRules(
            const RegularizationRules& base_rules,
            Real force_scale = 1.0_r,
            Real min_scale = 1.0_r
        );

        explicit ScaledRegularizationRules(
            PerConfigTypeSigmas base_sigmas,
            Real force_scale = 1.0_r,
            Real min_scale = 1.0_r
        );

        explicit ScaledRegularizationRules(
            Real energy_sigma_per_atom = 0.001,
            Real force_component_sigma = 0.05,
            Real virials_iso_sigma_per_atom = 0.1,
            Real virials_aniso_sigmas_per_atom = 0.02,
            Real force_scale = 1.0_r,
            Real min_scale = 1.0_r
        );

        Regularization determine(const Atoms& atoms) const override;

        ScaledRegularizationRules* clone() const override {
            return new ScaledRegularizationRules(*this);
        }

        Real getForceScale() const { return force_scale; }
        Real getMinScale() const { return min_scale; }
        const RegularizationRules* getBaseRules() const { return base_rules.get(); }

    private:
        std::shared_ptr<RegularizationRules> base_rules;
        Real force_scale;
        Real min_scale;
    };

}

#endif
