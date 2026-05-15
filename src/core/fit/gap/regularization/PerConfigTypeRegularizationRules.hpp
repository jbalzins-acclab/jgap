#ifndef JGAP_PERCONFIGTYPEREGULARIZATIONRULES_HPP
#define JGAP_PERCONFIGTYPEREGULARIZATIONRULES_HPP

#include "RegularizationRules.hpp"
#include "core/atomic/energy/Virials.hpp"
#include "core/atomic/geometry/Vector3.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class PerConfigTypeRegularizationRules : public RegularizationRules {
    public:
        PerConfigTypeRegularizationRules(Real energy_sigma_per_atom,
                                         Real force_component_sigma,
                                         Real virials_iso_sigma_per_atom,
                                         Real virials_aniso_sigmas_per_atom,
                                         const std::map<std::string, Real> &exact,
                                         const std::map<std::string, Real> &contains
        );

        void fillSigmas(Regularization &sigmas, const Atoms &atoms, const NeighbourList &neighbour_list) override;

    private:
        Real default_energy_per_atom;
        Vector3 default_force;
        Virials default_virials_per_atom;

        std::map<std::string, Real> exact;
        std::map<std::string, Real> contains; // order-sensitive
    };
}

#endif