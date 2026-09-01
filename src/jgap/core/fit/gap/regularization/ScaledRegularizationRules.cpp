#include "ScaledRegularizationRules.hpp"
#include <algorithm>
#include "PerConfigTypeRegularizationRules.hpp"
#include "SimpleRegularizationRules.hpp"

namespace jgap {

    ScaledRegularizationRules::ScaledRegularizationRules(
        std::shared_ptr<RegularizationRules> base_rules,
        const Real force_scale,
        const Real min_scale
    ) :
        base_rules(std::move(base_rules)),
        force_scale(force_scale),
        min_scale(min_scale) {}

    ScaledRegularizationRules::ScaledRegularizationRules(
        const RegularizationRules& base_rules,
        const Real force_scale,
        const Real min_scale
    ) :
        base_rules(base_rules.clone()),
        force_scale(force_scale),
        min_scale(min_scale) {}

    ScaledRegularizationRules::ScaledRegularizationRules(
        PerConfigTypeSigmas base_sigmas,
        const Real force_scale,
        const Real min_scale
    ) :
        base_rules(std::make_shared<PerConfigTypeRegularizationRules>(std::move(base_sigmas))),
        force_scale(force_scale),
        min_scale(min_scale) {}

    ScaledRegularizationRules::ScaledRegularizationRules(
        const Real energy_sigma_per_atom,
        const Real force_component_sigma,
        const Real virials_iso_sigma_per_atom,
        const Real virials_aniso_sigmas_per_atom,
        const Real force_scale,
        const Real min_scale
    ) :
        base_rules(std::make_shared<SimpleRegularizationRules>(
            energy_sigma_per_atom,
            force_component_sigma,
            virials_iso_sigma_per_atom,
            virials_aniso_sigmas_per_atom
        )),
        force_scale(force_scale),
        min_scale(min_scale) {}

    Regularization ScaledRegularizationRules::determine(const Atoms& atoms) const {
        Regularization sigmas = base_rules ? base_rules->determine(atoms) : Regularization{};

        Real max_force = 0.0_r;
        if (const auto& forces = atoms.getForces(); forces.has_value() && !forces->empty()) {
            for (const auto& f: *forces) {
                max_force = std::max(max_force, f.norm());
            }
        }

        const Real scale = std::max(min_scale, 1.0_r + force_scale * max_force);

        if (sigmas.energy.has_value()) {
            *sigmas.energy *= scale;
        }

        if (sigmas.forces.has_value()) {
            for (auto& f: *sigmas.forces) {
                f *= scale;
            }
        }

        if (sigmas.virials.has_value()) {
            *sigmas.virials = *sigmas.virials * scale;
        }

        return sigmas;
    }

} // namespace jgap
