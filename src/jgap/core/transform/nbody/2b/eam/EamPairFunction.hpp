#ifndef JGAP_EAMPAIRFUNCTION_HPP
#define JGAP_EAMPAIRFUNCTION_HPP

#include <map>
#include <set>
#include "../TwoBodyTransformation.hpp"
#include "jgap/core/ValuePtr.hpp"
#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/core/atomic/neighbours/NeighbourLists.hpp"
#include "jgap/core/potentials/gap/component/ManyBodyGapComponent.hpp"

namespace jgap {
    enum class EamMode { FSsym, FSgen, EAM, Blind };

    class EamPairFunction : public TwoBodyTransformation<1> {
    public:
        static std::map<Species, ValuePtr<NBodyAggregator<1>>> createAggregators(
            const ValuePtr<TwoBodyTransformation<1>>& base_pf, const std::vector<Atoms>& training_data, EamMode mode);

        template<CKernelOfDim<1> TKernel>
        static std::vector<ManyBodyGapComponent<1, TKernel>> createComponents(
            const ValuePtr<TwoBodyTransformation<1>>& base_pf, const TKernel& kernel, const Sparsifier<1>& sparsifier,
            const std::vector<Atoms>& training_data, EamMode mode, const std::vector<Real>& optional_coeffs = {});

        Cutoffs getCutoffs() const override { return Cutoffs{{2, cutoff}}; }

        void setPrefactor(Real p) { prefactor = p; }

        Real getPrefactor() const { return prefactor; }

        Real getCutoff() const { return cutoff; }

        EamPairFunction* clone() const override = 0;

    protected:
        Real cutoff;
        Real prefactor;

        EamPairFunction(Real cutoff = 0.0, Real prefactor = 1.0) : cutoff(cutoff), prefactor(prefactor) {}
    };

    template<CKernelOfDim<1> TKernel>
    std::vector<ManyBodyGapComponent<1, TKernel>> EamPairFunction::createComponents(
        const ValuePtr<TwoBodyTransformation<1>>& base_pf, const TKernel& kernel, const Sparsifier<1>& sparsifier,
        const std::vector<Atoms>& training_data, EamMode mode, const std::vector<Real>& optional_coeffs) {
        auto aggregators = createAggregators(base_pf, training_data, mode);
        std::vector<ManyBodyGapComponent<1, TKernel>> components;

        for (auto& [central_species, aggregator]: aggregators) {
            components.emplace_back(std::move(aggregator), kernel, sparsifier, training_data, optional_coeffs);
        }

        return components;
    }
}

#endif
