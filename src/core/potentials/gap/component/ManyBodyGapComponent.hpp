#ifndef JGAP_MANYBODYGAPCOMPONENT_HPP
#define JGAP_MANYBODYGAPCOMPONENT_HPP

#include "GapComponent.hpp"
#include "core/kernels/Kernel.hpp"

namespace jgap {

    template<typename TKernel> requires CKernel<TKernel>
    class ManyBodyGapComponent : public GapComponent {
    public:
        static constexpr size_t Dim = TKernel::Dim;
        static constexpr size_t ClusterSize = TTransformer::ClusterSize;

        ManyBodyGapComponent(Species central_atom_species) {}

        AtomicQuantities covariate(const NeighbourList &neighbour_list) const override {



        }

        void blabla(const NeighbourList &neighbour_list) const {
            auto clusters = neighbour_list.findAllClusters<Dim>(some_species_set);


            for (auto cluster: clusters) {
                DescriptorAndDerivatives<Dim, ClusterSize> contribution
                    = transformation.evaluateAndDifferentiate(cluster);


            }


            for (size_t sparse_idx = 0; sparse_idx < sparse_points.size(); sparse_idx++) {
                const auto [K, gradK_wrt_q] = kernel.valueAndGradient(
                    sparse_points[sparse_idx].value, descriptor.value
                    );
            }
        }

        Matrix sparseToSparseCovariance() const override {
            Matrix result(sparse_points.size(), sparse_points.size());

            for (size_t i = 0; i < sparse_points.size(); i++) {
                for (size_t j = i; j < sparse_points.size(); j++) {
                    result(i, j) = kernel.value(sparse_points[i].value, sparse_points[j].value);
                    result(j, i) = result(i, j);
                }
            }

            return result;
        }

        size_t nSparsePoints() const override {
            return sparse_points.size();
        }

        Cutoffs getCutoffs() const override {

        }

    private:
        TKernel kernel;

        std::vector<Descriptor<Dim>> sparse_points;
        SpeciesSet central_atom_species;
    };
}

#endif
