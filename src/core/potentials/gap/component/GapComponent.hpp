#ifndef JGAP_GAPCOMPONENT_HPP
#define JGAP_GAPCOMPONENT_HPP

#include <vector>

#include "core/Matrix.hpp"
#include "core/atomic/Atoms.hpp"
#include "core/atomic/energy/AtomicQuantities.hpp"
#include "core/atomic/energy/AtomicQuantity.hpp"
#include "core/atomic/neighbours/NeighbourList.hpp"
#include "core/transform/ClusterTransformation.hpp"
#include "io/Serializable.hpp"

namespace jgap {

    class GapComponent {
    public:
        using Ptr = std::unique_ptr<GapComponent>;

        virtual ~GapComponent() = default;
        virtual std::optional<AtomicQuantities> covariate(const NeighbourList& neighbour_list) const = 0;
        virtual Matrix sparseToSparseCovariance() const = 0;
        virtual size_t nSparsePoints() const = 0;
        virtual Cutoffs getCutoffs() const = 0;

        template <std::forward_iterator It>
        void setCoefficients(It& iter);

        void setCoefficients(const std::vector<Real>& new_coeff);

        Real getCutoff() const { return getCutoffs().maxOverall(); }

        AtomicQuantity energy(const Atoms& atoms) const;
        AtomicQuantity energy(const NeighbourList& neighbour_list) const;

    private:
        std::vector<Real> coefficients{};
    };

    template<std::forward_iterator It>
    void GapComponent::setCoefficients(It &iter) {
        coefficients.resize(nSparsePoints());
        for (size_t i{}; i < nSparsePoints(); i++, iter++) {
            coefficients[i] = *iter;
        }
    }

    /*
    template<size_t Dim, size_t Dependencies>
    class GapComponent : public GapComponent {
    public:
        //static constexpr size_t DescriptorDim = 4;
        //static constexpr size_t DescriptorDependencies = 3;
        using TTransformerGroup = TransformerGroup<Dim, Dependencies>;
        using TKernelGroup = SparseGroup<Dim, Dependencies>;

        GapComponent(TTransformerGroup transformation, TKernelGroup kernel_group)
            : transformation(transformation), kernel_group(kernel_group) {
        }

        size_t nSparsePoints() const override {
            return kernel_group.nSparsePoints();
        }

        std::vector<Matrix> sparseToSparseCovariance() const override {
            return kernel_group.sparseToSparseCovariance();
        }

        Cutoffs getCutoffs() const override {
            return transformation.getCutoffs();
        }

        std::vector<AtomicQuantity> covariate(const NeighbourList& neighbour_list) const override {
            auto descriptors = transformation.template transform<CalculationType::WithGradients>(neighbour_list);
            return kernel_group->covariate(neighbour_list.nAtoms(), descriptors);
        }

    private:
        TTransformerGroup transformation;
        std::shared_ptr<TKernelGroup> kernel_group;
    };
    */
}

#endif
