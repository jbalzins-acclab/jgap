#ifndef JGAP_NBODYGAPCOMPONENT_HPP
#define JGAP_NBODYGAPCOMPONENT_HPP

#include "core/atomic/energy/AtomicQuantities.hpp"
#include "core/kernels/Kernel.hpp"
#include "core/transform/ClusterTransformation.hpp"
#include "GapComponent.hpp"
#include "core/Matrix.hpp"

namespace jgap {

    template<size_t Dim, size_t ClusterSize, ClusterTypes ClusterType>
    class NBodyGapComponent : public GapComponent {
    public:
        static constexpr size_t Dependencies = Cluster<ClusterSize>::NSeparations;

        NBodyGapComponent(const SpeciesSet<ClusterSize, ClusterType> species,
                          std::unique_ptr<ClusterTransformation<Dim, ClusterSize> > transformation,
                          std::unique_ptr<Kernel<Dim> > kernel,
                          std::vector<Descriptor<Dim> > sparse_points,
                          const std::vector<Real> &optional_coeffs = {})
            : species(species),
              symmetry_factor(transformation->symmetryFactor()),
              transformation(std::move(transformation)),
              kernel(std::move(kernel)),
              sparse_points(std::move(sparse_points)) {

            if (!optional_coeffs.empty()) {
                setCoefficients(optional_coeffs);
            }
        }

        std::optional<AtomicQuantities> covariate(const NeighbourList& neighbour_list) const override {

            auto clusters = neighbour_list.findAllClusters<ClusterSize, ClusterType>(species);
            if (clusters.empty()) {
                return std::nullopt;
            }

            AtomicQuantities result(nSparsePoints(), neighbour_list.nAtoms());

            for (const auto& cluster: clusters) {
                auto descriptor = transformation->evaluateAndDifferentiate(cluster);

                for (size_t sparse_idx = 0; sparse_idx < sparse_points.size(); sparse_idx++) {
                    auto [K, gradK_wrt_q] = kernel->valueAndGradient(
                        sparse_points[sparse_idx].value, descriptor.value
                    );

                    K *= symmetry_factor;
                    for(auto& grad_val : gradK_wrt_q) {
                        grad_val *= symmetry_factor;
                    }

                    result.energy(sparse_idx) += K;

                    for (size_t i = 0; i < ClusterSize; i++) {
                        for (size_t j = i + 1; j < ClusterSize; j++) {
                            const auto sep_idx = flattened_index(i, j);
                            const auto& separation = cluster.between(i, j);
                            const auto& derivatives_wrt_r_norms = descriptor.derivatives[sep_idx];

                            Real dK_drij_norm = 0.0;
                            for (size_t dim = 0; dim < Dim; dim++) {
                                dK_drij_norm += derivatives_wrt_r_norms[dim] * gradK_wrt_q[dim];
                            }

                            result.force(sparse_idx, cluster.atom_indexes[i]) += separation.direction * dK_drij_norm;
                            result.force(sparse_idx, cluster.atom_indexes[j]) -= separation.direction * dK_drij_norm;

                            result.virials(sparse_idx) += separation.virials * dK_drij_norm;
                        }
                    }
                }
            }

            return result;
        }

        Matrix sparseToSparseCovariance() const override {
            Matrix result(nSparsePoints(), nSparsePoints());
            for (size_t i = 0; i < nSparsePoints(); i++) {
                for (size_t j = i; j < nSparsePoints(); j++) {
                    result(i, j) = kernel->value(sparse_points[i].value, sparse_points[j].value);
                    result(j, i) = result(i, j);
                }
            }
            return result;
        }

        size_t nSparsePoints() const override {
            return sparse_points.size();
        }

        Cutoffs getCutoffs() const override {
            return transformation->getCutoffs();
        }

    public:
        SpeciesSet<ClusterSize, ClusterType> species;
        Real symmetry_factor;
        std::unique_ptr<ClusterTransformation<Dim, ClusterSize>> transformation;
        std::unique_ptr<Kernel<Dim>> kernel;
        std::vector<Descriptor<Dim>> sparse_points;
    };
}

#endif