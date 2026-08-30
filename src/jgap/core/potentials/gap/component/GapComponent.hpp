#ifndef JGAP_GAPCOMPONENT_HPP
#define JGAP_GAPCOMPONENT_HPP

#include <set>
#include <vector>

#include "jgap/core/Matrix.hpp"
#include "jgap/core/ValuePtr.hpp"
#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/core/atomic/energy/AtomicQuantities.hpp"
#include "jgap/core/atomic/energy/AtomicQuantity.hpp"
#include "jgap/core/atomic/neighbours/NeighbourLists.hpp"
#include "jgap/core/atomic/species/Species.hpp"
#include "jgap/core/tabulation/TabulationData.hpp"

namespace jgap {

    class GapComponent {
    public:
        virtual ~GapComponent() = default;
        virtual std::optional<AtomicQuantities> covariate(const NeighbourLists& neighbour_list) const = 0;
        virtual Matrix<RowMajor> sparseToSparseCovariance() const = 0;
        virtual size_t nSparsePoints() const = 0;
        virtual Cutoffs getCutoffs() const = 0;
        virtual std::set<Species> nonZeroCovarianceFor() const = 0;

        virtual GapComponent* clone() const = 0;

        template<std::forward_iterator It>
        void setCoefficients(It& iter);

        void setCoefficients(const std::vector<Real>& new_coeff);
        const std::vector<Real>& getCoefficients() const { return coefficients; }

        Real getCutoff() const { return getCutoffs().maxOverall(); }

        AtomicQuantity energy(const Atoms& atoms) const;
        AtomicQuantity energy(const NeighbourLists& neighbour_list) const;

        virtual void tabulate(TabulationData& tables) const = 0;

    private:
        std::vector<Real> coefficients{};
    };

    template<std::forward_iterator It>
    void GapComponent::setCoefficients(It& iter) {
        coefficients.resize(nSparsePoints());
        for (size_t i{}; i < nSparsePoints(); i++, iter++) {
            coefficients[i] = *iter;
        }
    }

    static_assert(Cloneable<GapComponent>);
}

#endif
