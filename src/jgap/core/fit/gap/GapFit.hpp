#ifndef JGAP_GAPFIT_HPP
#define JGAP_GAPFIT_HPP

#include "jgap/core/potentials/gap/GapPotential.hpp"
#include "jgap/core/atomic/Atoms.hpp"
#include <vector>
#include <memory>

#include "regularization/Regularization.hpp"
#include "regularization/RegularizationRules.hpp"

namespace jgap {
    class GapFit {
    public:
        virtual ~GapFit() = default;

        void fit(GapPotential& to_be_fit,
                 const std::vector<Atoms>& training_data,
                 const RegularizationRules& regularization_rules,
                 const std::vector<Regularization>& sigmas = {});

    protected:
        struct EnergyData {
            std::optional<Real> energy;
            std::optional<Virials> virials;
            std::optional<std::vector<Vector3>> forces;
        };

        virtual std::vector<Real> findCoefficients(std::vector<ValuePtr<GapComponent>>& gap_components,
                                                   const std::vector<Atoms>& training_data,
                                                   std::vector<EnergyData>& energies_without_external,
                                                   std::vector<Regularization>& sigmas_inverse) = 0;
    };
}

#endif
