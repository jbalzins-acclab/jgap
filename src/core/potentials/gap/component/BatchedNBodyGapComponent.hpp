#ifndef JGAP_BATCHEDNBODYGAPCOMPONENT_HPP
#define JGAP_BATCHEDNBODYGAPCOMPONENT_HPP

#include <Eigen/Dense>
#include "GapComponent.hpp"
#include "core/kernels/SquaredExpKernel.hpp"
#include "core/transform/ClusterTransformation.hpp"
#include "core/atomic/energy/AtomicQuantities.hpp"

#ifdef JGAP_USE_MPS
#  ifdef __OBJC__
// Compiling as Objective-C++, import full headers.
#    import <Metal/Metal.h>
#    import <MetalPerformanceShaders/MetalPerformanceShaders.h>
#    import <Foundation/Foundation.h>
// Also define C-style typedefs for consistency.
typedef id<MTLDevice> id_MTLDevice;
typedef id<MTLCommandQueue> id_MTLCommandQueue;
typedef id<MTLComputePipelineState> id_MTLComputePipelineState;
typedef id<MTLBuffer> id_MTLBuffer;
#  else
// Compiling as C++, use opaque pointers.
typedef struct MTLDevice_struct* id_MTLDevice;
typedef struct MTLCommandQueue_struct* id_MTLCommandQueue;
typedef struct MTLComputePipelineState_struct* id_MTLComputePipelineState;
typedef struct MTLBuffer_struct* id_MTLBuffer;
#  endif
#endif

namespace jgap {

#ifdef JGAP_USE_MPS
namespace detail {
    struct MtlSharedBuf {
#ifdef __OBJC__
        id_MTLBuffer buf  = nil;
#else
        id_MTLBuffer buf  = nullptr;
#endif
        size_t       rows = 0;
        size_t       cols = 1;

        void alloc(id_MTLDevice dev, size_t r, size_t c = 1);
        float* data() const;
    };
}
#endif // JGAP_USE_MPS

template<size_t Dim, size_t ClusterSize, ClusterTypes ClusterType,
         size_t ExpDim, size_t CutDim>
requires (Dim == ExpDim + CutDim && CutDim == 1)
class BatchedNBodyGapComponent : public GapComponent {
public:
    static constexpr size_t NSep = Cluster<ClusterSize>::NSeparations;

    using EigenMat = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using EigenVec = Eigen::Matrix<Real, Eigen::Dynamic, 1>;

    BatchedNBodyGapComponent(
        SpeciesSet<ClusterSize, ClusterType> species,
        std::unique_ptr<ClusterTransformation<Dim, ClusterSize>> transformation,
        SquaredExpKernel<ExpDim, CutDim> kernel,
        std::vector<Descriptor<Dim>> sparse_points,
        const std::vector<Real>& optional_coeffs = {}
    );

    ~BatchedNBodyGapComponent();

    std::optional<AtomicQuantities> covariate(const NeighbourList& nl) const override;
    Matrix sparseToSparseCovariance() const override;
    size_t nSparsePoints() const override;
    Cutoffs getCutoffs() const override;

public:
    SpeciesSet<ClusterSize, ClusterType> species;
    Real symmetry_factor;
    std::unique_ptr<ClusterTransformation<Dim, ClusterSize>> transformation;
    SquaredExpKernel<ExpDim, CutDim> kernel;
    std::vector<Descriptor<Dim>> sparse_points;

private:
    static constexpr size_t CLUSTER_BATCH_SIZE =
#ifdef JGAP_USE_MPS
        512;
#else
        256;
#endif

    EigenMat QS_;
    EigenMat QS_se_;
    EigenVec norm_QS_;
    EigenVec cutoff_s_;

    void precomputeSparse();

    static constexpr std::array<std::pair<size_t, size_t>, NSep> computeAtomPairs() {
        std::array<std::pair<size_t, size_t>, NSep> pairs{};
        size_t idx = 0;
        for (size_t i = 0; i < ClusterSize; ++i)
            for (size_t j = i + 1; j < ClusterSize; ++j)
                pairs[idx++] = {i, j};
        return pairs;
    }

#ifndef JGAP_USE_MPS
    mutable EigenMat buf_Q_;
    mutable EigenMat buf_Q_se_;
    mutable EigenMat buf_K_;
    mutable EigenMat buf_K_ncs_;
    mutable EigenMat buf_fm_;
    mutable EigenMat buf_dQ_se_;
    mutable std::array<EigenMat, NSep> buf_dQ_;

    void allocateWorkBuffers();

    template<size_t NPairs>
    void covariateCpu(
        const auto& clusters, size_t N_c, size_t N_s,
        const auto& inv_l2, Real pf,
        const std::array<std::pair<size_t,size_t>, NPairs>& pairs,
        AtomicQuantities& result) const;
#endif

#ifdef JGAP_USE_MPS
#  ifdef __OBJC__
    id_MTLDevice               mtl_device_     = nil;
    id_MTLCommandQueue         mtl_queue_      = nil;
    id_MTLComputePipelineState ps_exp_scale_   = nil;
    id_MTLComputePipelineState ps_row_scale_   = nil;
    id_MTLComputePipelineState ps_col_scale_   = nil;
    id_MTLComputePipelineState ps_col_sub_     = nil;
    id_MTLComputePipelineState ps_fm_fma_      = nil;
#  else
    id_MTLDevice               mtl_device_     = nullptr;
    id_MTLCommandQueue         mtl_queue_      = nullptr;
    id_MTLComputePipelineState ps_exp_scale_   = nullptr;
    id_MTLComputePipelineState ps_row_scale_   = nullptr;
    id_MTLComputePipelineState ps_col_scale_   = nullptr;
    id_MTLComputePipelineState ps_col_sub_     = nullptr;
    id_MTLComputePipelineState ps_fm_fma_      = nullptr;
#  endif

    mutable detail::MtlSharedBuf gpu_Q_;
    mutable detail::MtlSharedBuf gpu_Q_se_;
    mutable detail::MtlSharedBuf gpu_norm_Q_;
    mutable detail::MtlSharedBuf gpu_cutoff_c_;
    mutable detail::MtlSharedBuf gpu_K_;
    mutable detail::MtlSharedBuf gpu_K_ncs_;
    mutable detail::MtlSharedBuf gpu_fm_;
    mutable detail::MtlSharedBuf gpu_dQ_se_;
    mutable detail::MtlSharedBuf gpu_self_;
    mutable detail::MtlSharedBuf gpu_dQ_cut_;
    mutable std::array<detail::MtlSharedBuf, NSep> gpu_dQ_;

    detail::MtlSharedBuf gpu_QS_se_;
    detail::MtlSharedBuf gpu_norm_QS_;
    detail::MtlSharedBuf gpu_cutoff_s_;

    void initMetal();
    void allocateGpuBuffers();

    template<size_t NPairs>
    void covariateGpu(
        const auto& clusters, size_t N_c, size_t N_s,
        const auto& inv_l2, Real pf,
        const std::array<std::pair<size_t,size_t>, NPairs>& pairs,
        AtomicQuantities& result) const;
#endif
};

// ============================================================================
// Implementation Details
// ============================================================================

template<size_t Dim, size_t ClusterSize, ClusterTypes ClusterType, size_t ExpDim, size_t CutDim>
requires (Dim == ExpDim + CutDim && CutDim == 1)
BatchedNBodyGapComponent<Dim, ClusterSize, ClusterType, ExpDim, CutDim>::BatchedNBodyGapComponent(
    SpeciesSet<ClusterSize, ClusterType> species,
    std::unique_ptr<ClusterTransformation<Dim, ClusterSize>> transformation,
    SquaredExpKernel<ExpDim, CutDim> kernel,
    std::vector<Descriptor<Dim>> sparse_points,
    const std::vector<Real>& optional_coeffs
)
    : species(std::move(species)),
      symmetry_factor(transformation->symmetryFactor()),
      transformation(std::move(transformation)),
      kernel(std::move(kernel)),
      sparse_points(std::move(sparse_points))
{
    if (!optional_coeffs.empty()) setCoefficients(optional_coeffs);
    precomputeSparse();
#ifdef JGAP_USE_MPS
    initMetal();
    allocateGpuBuffers();
#else
    allocateWorkBuffers();
#endif
}

template<size_t Dim, size_t ClusterSize, ClusterTypes ClusterType, size_t ExpDim, size_t CutDim>
requires (Dim == ExpDim + CutDim && CutDim == 1)
BatchedNBodyGapComponent<Dim, ClusterSize, ClusterType, ExpDim, CutDim>::~BatchedNBodyGapComponent() {
    // With ARC, we don't need to manually release Metal objects.
}

template<size_t Dim, size_t ClusterSize, ClusterTypes ClusterType, size_t ExpDim, size_t CutDim>
requires (Dim == ExpDim + CutDim && CutDim == 1)
std::optional<AtomicQuantities> BatchedNBodyGapComponent<Dim, ClusterSize, ClusterType, ExpDim, CutDim>::covariate(const NeighbourList& nl) const {
    auto clusters = nl.findAllClusters<ClusterSize, ClusterType>(species);
    if (clusters.empty()) return std::nullopt;

    const size_t N_c = clusters.size();
    const size_t N_s = nSparsePoints();
    AtomicQuantities result(N_s, nl.nAtoms());

    const auto& inv_l2 = kernel.getInverseLengthScalesSquared();
    const Real  pf     = kernel.getPrefactor();

    static constexpr auto pairs = computeAtomPairs();

#ifdef JGAP_USE_MPS
    covariateGpu(clusters, N_c, N_s, inv_l2, pf, pairs, result);
#else
    covariateCpu(clusters, N_c, N_s, inv_l2, pf, pairs, result);
#endif
    return result;
}

template<size_t Dim, size_t ClusterSize, ClusterTypes ClusterType, size_t ExpDim, size_t CutDim>
requires (Dim == ExpDim + CutDim && CutDim == 1)
Matrix BatchedNBodyGapComponent<Dim, ClusterSize, ClusterType, ExpDim, CutDim>::sparseToSparseCovariance() const {
    const size_t n = nSparsePoints();
    Matrix result(n, n);
    for (size_t i = 0; i < n; ++i)
        for (size_t j = i; j < n; ++j) {
            result(i, j) = kernel.value(sparse_points[i].value, sparse_points[j].value);
            result(j, i) = result(i, j);
        }
    return result;
}

template<size_t Dim, size_t ClusterSize, ClusterTypes ClusterType, size_t ExpDim, size_t CutDim>
requires (Dim == ExpDim + CutDim && CutDim == 1)
size_t BatchedNBodyGapComponent<Dim, ClusterSize, ClusterType, ExpDim, CutDim>::nSparsePoints() const {
    return sparse_points.size();
}

template<size_t Dim, size_t ClusterSize, ClusterTypes ClusterType, size_t ExpDim, size_t CutDim>
requires (Dim == ExpDim + CutDim && CutDim == 1)
Cutoffs BatchedNBodyGapComponent<Dim, ClusterSize, ClusterType, ExpDim, CutDim>::getCutoffs() const {
    return transformation->getCutoffs();
}

template<size_t Dim, size_t ClusterSize, ClusterTypes ClusterType, size_t ExpDim, size_t CutDim>
requires (Dim == ExpDim + CutDim && CutDim == 1)
void BatchedNBodyGapComponent<Dim, ClusterSize, ClusterType, ExpDim, CutDim>::precomputeSparse() {
    const size_t n      = sparse_points.size();
    const auto&  inv_l2 = kernel.getInverseLengthScalesSquared();

    QS_.resize(n, Dim);
    for (size_t si = 0; si < n; ++si)
        for (size_t d = 0; d < Dim; ++d)
            QS_(si, d) = sparse_points[si].value[d];

    QS_se_ = QS_.leftCols(ExpDim);
    for (size_t d = 0; d < ExpDim; ++d)
        QS_se_.col(d) *= std::sqrt(inv_l2[d]);

    norm_QS_  = QS_se_.rowwise().squaredNorm();
    cutoff_s_ = QS_.col(ExpDim);
}

#ifndef JGAP_USE_MPS
template<size_t Dim, size_t ClusterSize, ClusterTypes ClusterType, size_t ExpDim, size_t CutDim>
requires (Dim == ExpDim + CutDim && CutDim == 1)
void BatchedNBodyGapComponent<Dim, ClusterSize, ClusterType, ExpDim, CutDim>::allocateWorkBuffers() {
    const size_t B  = CLUSTER_BATCH_SIZE;
    const size_t Ns = sparse_points.size();
    buf_Q_.resize(B, Dim);
    buf_Q_se_.resize(B, ExpDim);
    buf_K_.resize(B, Ns);
    buf_K_ncs_.resize(B, Ns);
    buf_fm_.resize(B, Ns);
    buf_dQ_se_.resize(B, ExpDim);
    for (auto& m : buf_dQ_) m.resize(B, Dim);
}

template<size_t Dim, size_t ClusterSize, ClusterTypes ClusterType, size_t ExpDim, size_t CutDim>
requires (Dim == ExpDim + CutDim && CutDim == 1)
template<size_t NPairs>
void BatchedNBodyGapComponent<Dim, ClusterSize, ClusterType, ExpDim, CutDim>::covariateCpu(
    const auto& clusters, size_t N_c, size_t N_s,
    const auto& inv_l2, Real pf,
    const std::array<std::pair<size_t,size_t>, NPairs>& pairs,
    AtomicQuantities& result) const
{
    for (size_t batch_start = 0; batch_start < N_c; batch_start += CLUSTER_BATCH_SIZE) {
        const size_t B = std::min(CLUSTER_BATCH_SIZE, N_c - batch_start);

        for (size_t ci = 0; ci < B; ++ci) {
            auto desc = transformation->evaluateAndDifferentiate(clusters[batch_start + ci]);
            for (size_t d = 0; d < Dim; ++d)
                buf_Q_(ci, d) = desc.value[d];
            for (size_t sep = 0; sep < NSep; ++sep)
                for (size_t d = 0; d < Dim; ++d)
                    buf_dQ_[sep](ci, d) = desc.derivatives[sep][d];
        }

        auto Q        = buf_Q_.topRows(B);
        auto Q_se_buf = buf_Q_se_.topRows(B);
        Q_se_buf = Q.leftCols(ExpDim);
        for (size_t d = 0; d < ExpDim; ++d)
            Q_se_buf.col(d) *= std::sqrt(inv_l2[d]);

        const EigenVec norm_Q = Q_se_buf.rowwise().squaredNorm();

        auto dist2 = buf_K_.topRows(B);
        dist2.noalias() = -2.0 * (Q_se_buf * QS_se_.transpose());
        dist2.colwise() += norm_Q;
        dist2.rowwise() += norm_QS_.transpose();
        dist2 = dist2.cwiseMax(0.0);
        dist2.array() = (dist2.array() * (-0.5)).exp() * pf;

        const EigenVec cutoff_c = Q.col(ExpDim);
        dist2.array().rowwise() *= (cutoff_s_ * symmetry_factor).transpose().array();
        buf_K_ncs_.topRows(B) = dist2;
        dist2.array().colwise() *= cutoff_c.array();
        const auto& K            = dist2;
        const auto  K_no_cut_sym = buf_K_ncs_.topRows(B);

        const EigenVec energy_vec = K.colwise().sum();
        for (size_t si = 0; si < N_s; ++si)
            result.energy(si) += energy_vec(si);

        for (size_t sep = 0; sep < NSep; ++sep) {
            const auto [ai, aj] = pairs[sep];
            auto dQ_sep = buf_dQ_[sep].topRows(B);
            auto dQ_se  = buf_dQ_se_.topRows(B);
            dQ_se = dQ_sep.leftCols(ExpDim);
            for (size_t d = 0; d < ExpDim; ++d)
                dQ_se.col(d) *= inv_l2[d];

            const EigenVec self =
                (dQ_se.array() * Q_se_buf.array()).rowwise().sum();

            auto fm = buf_fm_.topRows(B);
            fm.noalias() = K.cwiseProduct(
                (dQ_se * QS_se_.transpose()).colwise() - self);
            fm.array() +=
                K_no_cut_sym.array().colwise() * dQ_sep.col(ExpDim).array();

            for (size_t ci = 0; ci < B; ++ci) {
                const size_t atom_i  = clusters[batch_start + ci].atom_indexes[ai];
                const size_t atom_j  = clusters[batch_start + ci].atom_indexes[aj];
                const auto&  sep_geo = clusters[batch_start + ci].separations[sep];
                for (size_t si = 0; si < N_s; ++si) {
                    const Real fmv = fm(ci, si);
                    result.force(si, atom_i) += sep_geo.direction * fmv;
                    result.force(si, atom_j) -= sep_geo.direction * fmv;
                    result.virials(si)        += sep_geo.virials   * fmv;
                }
            }
        }
    }
}
#endif // !JGAP_USE_MPS

} // namespace jgap

#endif // JGAP_BATCHEDNBODYGAPCOMPONENT_HPP