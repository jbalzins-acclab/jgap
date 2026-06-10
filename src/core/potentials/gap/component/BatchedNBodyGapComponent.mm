#include "BatchedNBodyGapComponent.hpp"

#ifdef JGAP_USE_MPS

namespace jgap {
namespace detail {

// ─── Custom Metal kernel source ───────────────────────────────────────────────
static constexpr const char* kKernelSrc = R"MSL(
#include <metal_stdlib>
using namespace metal;

kernel void kernel_exp_scale(
    device float* A      [[buffer(0)]],
    constant float& scale[[buffer(1)]],
    uint id [[thread_position_in_grid]])
{
    A[id] = exp(A[id] * scale);
}

kernel void kernel_row_scale(
    device float*        A    [[buffer(0)]],
    constant float*      v    [[buffer(1)]],
    constant uint&       cols [[buffer(2)]],
    uint2 gid [[thread_position_in_grid]])
{
    A[gid.y * cols + gid.x] *= v[gid.x];
}

kernel void kernel_col_scale(
    device float*        A    [[buffer(0)]],
    constant float*      v    [[buffer(1)]],
    constant uint&       cols [[buffer(2)]],
    uint2 gid [[thread_position_in_grid]])
{
    A[gid.y * cols + gid.x] *= v[gid.y];
}

kernel void kernel_col_sub(
    device float*        A    [[buffer(0)]],
    constant float*      v    [[buffer(1)]],
    constant uint&       cols [[buffer(2)]],
    uint2 gid [[thread_position_in_grid]])
{
    A[gid.y * cols + gid.x] -= v[gid.y];
}

kernel void kernel_fm_fma(
    device       float* fm   [[buffer(0)]],
    constant     float* K    [[buffer(1)]],
    constant     float* ncs  [[buffer(2)]],
    constant     float* dcut [[buffer(3)]],
    constant     uint&  cols [[buffer(4)]],
    uint2 gid [[thread_position_in_grid]])
{
    const uint idx = gid.y * cols + gid.x;
    fm[idx] = K[idx] * fm[idx] + ncs[idx] * dcut[gid.y];
}
)MSL";

float* MtlSharedBuf::data() const { return static_cast<float*>([buf contents]); }

void MtlSharedBuf::alloc(id_MTLDevice dev, size_t r, size_t c) {
    rows = r; cols = c;
    const size_t bytes = r * c * sizeof(float);
    buf = [dev newBufferWithLength: bytes
                           options: MTLResourceStorageModeShared];
}

inline MPSMatrixDescriptor* makeDesc(NSUInteger rows, NSUInteger cols) {
    return [MPSMatrixDescriptor
        matrixDescriptorWithRows: rows
                         columns: cols
                        rowBytes: cols * sizeof(float)
                        dataType: MPSDataTypeFloat32];
}

inline MPSMatrix* asMPS(const MtlSharedBuf& b) {
    return [[MPSMatrix alloc]
        initWithBuffer: b.buf
            descriptor: makeDesc(b.rows, b.cols)];
}

void encodeElemwise1D(
    id<MTLComputeCommandEncoder> enc,
    id<MTLComputePipelineState>  pso,
    id<MTLBuffer>                buf,
    size_t                       nelem,
    const void*                  param,
    size_t                       paramBytes)
{
    [enc setComputePipelineState: pso];
    [enc setBuffer: buf offset: 0 atIndex: 0];
    [enc setBytes: param length: paramBytes atIndex: 1];
    const NSUInteger tpg = std::min((NSUInteger)nelem,
        pso.maxTotalThreadsPerThreadgroup);
    [enc dispatchThreads: MTLSizeMake(nelem, 1, 1)
        threadsPerThreadgroup: MTLSizeMake(tpg, 1, 1)];
}

void encodeElemwise2D(
    id<MTLComputeCommandEncoder> enc,
    id<MTLComputePipelineState>  pso,
    id<MTLBuffer>                A,
    id<MTLBuffer>                v,
    uint32_t                     rows,
    uint32_t                     cols)
{
    [enc setComputePipelineState: pso];
    [enc setBuffer: A offset: 0 atIndex: 0];
    [enc setBuffer: v offset: 0 atIndex: 1];
    [enc setBytes: &cols length: sizeof(uint32_t) atIndex: 2];
    const MTLSize grid = MTLSizeMake(cols, rows, 1);
    const MTLSize tpg  = MTLSizeMake(
        std::min((NSUInteger)16, (NSUInteger)cols),
        std::min((NSUInteger)16, (NSUInteger)rows), 1);
    [enc dispatchThreads: grid threadsPerThreadgroup: tpg];
}

} // namespace detail


template<size_t Dim, size_t ClusterSize, ClusterTypes ClusterType, size_t ExpDim, size_t CutDim>
requires (Dim == ExpDim + CutDim && CutDim == 1)
void BatchedNBodyGapComponent<Dim, ClusterSize, ClusterType, ExpDim, CutDim>::initMetal() {
    mtl_device_ = MTLCreateSystemDefaultDevice();
    if (!mtl_device_)
        throw std::runtime_error("No Metal device found — is this Apple Silicon?");
    mtl_queue_  = [mtl_device_ newCommandQueue];

    NSError* err = nil;
    id<MTLLibrary> lib = [mtl_device_
        newLibraryWithSource: [NSString stringWithUTF8String: detail::kKernelSrc]
                     options: nil
                       error: &err];
    if (!lib)
        throw std::runtime_error(
            std::string("Metal kernel compile failed: ") +
            [[err localizedDescription] UTF8String]);

    auto makePSO = [&](const char* name) {
        id<MTLFunction> fn = [lib newFunctionWithName:
            [NSString stringWithUTF8String: name]];
        id<MTLComputePipelineState> pso =
            [mtl_device_ newComputePipelineStateWithFunction: fn error: &err];
        if (!pso)
            throw std::runtime_error(
                std::string("PSO creation failed for ") + name + ": " +
                [[err localizedDescription] UTF8String]);
        return pso;
    };

    ps_exp_scale_ = makePSO("kernel_exp_scale");
    ps_row_scale_ = makePSO("kernel_row_scale");
    ps_col_scale_ = makePSO("kernel_col_scale");
    ps_col_sub_   = makePSO("kernel_col_sub");
    ps_fm_fma_    = makePSO("kernel_fm_fma");
}

template<size_t Dim, size_t ClusterSize, ClusterTypes ClusterType, size_t ExpDim, size_t CutDim>
requires (Dim == ExpDim + CutDim && CutDim == 1)
void BatchedNBodyGapComponent<Dim, ClusterSize, ClusterType, ExpDim, CutDim>::allocateGpuBuffers() {
    const size_t B  = CLUSTER_BATCH_SIZE;
    const size_t Ns = sparse_points.size();

    gpu_Q_.alloc(mtl_device_, B, Dim);
    gpu_Q_se_.alloc(mtl_device_, B, ExpDim);
    gpu_norm_Q_.alloc(mtl_device_, B);
    gpu_cutoff_c_.alloc(mtl_device_, B);
    gpu_K_.alloc(mtl_device_, B, Ns);
    gpu_K_ncs_.alloc(mtl_device_, B, Ns);
    gpu_fm_.alloc(mtl_device_, B, Ns);
    gpu_dQ_se_.alloc(mtl_device_, B, ExpDim);
    gpu_self_.alloc(mtl_device_, B);
    gpu_dQ_cut_.alloc(mtl_device_, B);
    for (auto& b : gpu_dQ_) b.alloc(mtl_device_, B, Dim);

    const auto& inv_l2 = kernel.getInverseLengthScalesSquared();
    gpu_QS_se_.alloc(mtl_device_, Ns, ExpDim);
    gpu_norm_QS_.alloc(mtl_device_, Ns);
    gpu_cutoff_s_.alloc(mtl_device_, Ns);

    Eigen::MatrixXf qsse = QS_se_.template cast<float>();
    std::memcpy(gpu_QS_se_.data(), qsse.data(), Ns * ExpDim * sizeof(float));

    Eigen::VectorXf nqs = norm_QS_.template cast<float>();
    std::memcpy(gpu_norm_QS_.data(), nqs.data(), Ns * sizeof(float));

    Eigen::VectorXf cs = cutoff_s_.template cast<float>();
    std::memcpy(gpu_cutoff_s_.data(), cs.data(), Ns * sizeof(float));
}

template<size_t Dim, size_t ClusterSize, ClusterTypes ClusterType, size_t ExpDim, size_t CutDim>
requires (Dim == ExpDim + CutDim && CutDim == 1)
template<size_t NPairs>
void BatchedNBodyGapComponent<Dim, ClusterSize, ClusterType, ExpDim, CutDim>::covariateGpu(
    const auto& clusters, size_t N_c, size_t N_s,
    const auto& inv_l2, Real pf,
    const std::array<std::pair<size_t,size_t>, NPairs>& pairs,
    AtomicQuantities& result) const
{
    const float pf_f   = static_cast<float>(pf);
    const float half_neg = -0.5f;

    std::vector<float> cutoff_s_sym(N_s);
    for (size_t si = 0; si < N_s; ++si)
        cutoff_s_sym[si] = static_cast<float>(cutoff_s_[si]) * symmetry_factor;

    for (size_t batch_start = 0; batch_start < N_c; batch_start += CLUSTER_BATCH_SIZE) {
        const size_t B = std::min(CLUSTER_BATCH_SIZE, N_c - batch_start);

        float* q_ptr    = gpu_Q_.data();
        float* norm_ptr = gpu_norm_Q_.data();
        float* cut_ptr  = gpu_cutoff_c_.data();

        for (size_t ci = 0; ci < B; ++ci) {
            auto desc = transformation->evaluateAndDifferentiate(
                clusters[batch_start + ci]);

            for (size_t d = 0; d < Dim; ++d)
                q_ptr[ci * Dim + d] = static_cast<float>(desc.value[d]);

            float norm_q = 0.0f;
            float* q_se_row = gpu_Q_se_.data() + ci * ExpDim;
            for (size_t d = 0; d < ExpDim; ++d) {
                const float v = static_cast<float>(desc.value[d]) *
                                std::sqrt(static_cast<float>(inv_l2[d]));
                q_se_row[d] = v;
                norm_q     += v * v;
            }
            norm_ptr[ci] = norm_q;
            cut_ptr[ci]  = static_cast<float>(desc.value[ExpDim]);

            for (size_t sep = 0; sep < NSep; ++sep) {
                float* dq_row = gpu_dQ_[sep].data() + ci * Dim;
                for (size_t d = 0; d < Dim; ++d)
                    dq_row[d] = static_cast<float>(desc.derivatives[sep][d]);
            }
        }

        id<MTLCommandBuffer> cmd = [mtl_queue_ commandBuffer];

        {
            auto* descA  = detail::makeDesc(B,   ExpDim);
            auto* descB  = detail::makeDesc(N_s, ExpDim);
            auto* descC  = detail::makeDesc(B,   N_s);

            auto* matA = detail::asMPS(gpu_Q_se_);
            auto* matB = detail::asMPS(gpu_QS_se_);
            auto* matC = detail::asMPS(gpu_K_);

            auto* gemm = [[MPSMatrixMultiplication alloc]
                initWithDevice:     mtl_device_
                transposeLeft:      NO
                transposeRight:     YES
                resultRows:         B
                resultColumns:      N_s
                interiorColumns:    ExpDim
                alpha:              -2.0
                beta:               0.0];

            [gemm encodeToCommandBuffer: cmd
                           leftMatrix:  matA
                          rightMatrix:  matB
                         resultMatrix:  matC];
        }

        [cmd commit];
        [cmd waitUntilCompleted];

        {
            float* k = gpu_K_.data();
            const float* nq  = gpu_norm_Q_.data();
            const float* nqs = gpu_norm_QS_.data();
            for (size_t ci = 0; ci < B; ++ci)
                for (size_t si = 0; si < N_s; ++si) {
                    float v = k[ci * N_s + si] + nq[ci] + nqs[si];
                    k[ci * N_s + si] = v < 0.0f ? 0.0f : v;
                }
        }

        cmd = [mtl_queue_ commandBuffer];
        {
            id<MTLComputeCommandEncoder> enc =
                [cmd computeCommandEncoder];

            const float neg_half = -0.5f;
            detail::encodeElemwise1D(enc, ps_exp_scale_, gpu_K_.buf,
                             B * N_s, &neg_half, sizeof(float));

            std::vector<float> cs_pf(N_s);
            for (size_t si = 0; si < N_s; ++si)
                cs_pf[si] = cutoff_s_sym[si] * pf_f;

            std::memcpy(gpu_cutoff_s_.data(), cs_pf.data(), N_s * sizeof(float));

            detail::encodeElemwise2D(enc, ps_row_scale_,
                             gpu_K_.buf, gpu_cutoff_s_.buf,
                             static_cast<uint32_t>(B),
                             static_cast<uint32_t>(N_s));

            [enc endEncoding];

            id<MTLBlitCommandEncoder> blit = [cmd blitCommandEncoder];
            [blit copyFromBuffer: gpu_K_.buf
                    sourceOffset: 0
                        toBuffer: gpu_K_ncs_.buf
               destinationOffset: 0
                            size: B * N_s * sizeof(float)];
            [blit endEncoding];

            id<MTLComputeCommandEncoder> enc2 = [cmd computeCommandEncoder];
            detail::encodeElemwise2D(enc2, ps_col_scale_,
                             gpu_K_.buf, gpu_cutoff_c_.buf,
                             static_cast<uint32_t>(B),
                             static_cast<uint32_t>(N_s));
            [enc2 endEncoding];
        }

        [cmd commit];
        [cmd waitUntilCompleted];

        {
            const float* k = gpu_K_.data();
            for (size_t si = 0; si < N_s; ++si) {
                float e = 0.0f;
                for (size_t ci = 0; ci < B; ++ci)
                    e += k[ci * N_s + si];
                result.energy(si) += static_cast<Real>(e);
            }
        }

        for (size_t sep = 0; sep < NSep; ++sep) {
            const auto [ai, aj] = pairs[sep];

            {
                float* dq_se = gpu_dQ_se_.data();
                float* self  = gpu_self_.data();
                float* dq_cut = gpu_dQ_cut_.data();

                const float* q_se = gpu_Q_se_.data();

                for (size_t ci = 0; ci < B; ++ci) {
                    float s = 0.0f;
                    for (size_t d = 0; d < ExpDim; ++d) {
                        const float v = static_cast<float>(
                            gpu_dQ_[sep].data()[ci * Dim + d]) *
                            static_cast<float>(inv_l2[d]);
                        dq_se[ci * ExpDim + d] = v;
                        s += v * q_se[ci * ExpDim + d];
                    }
                    self[ci]   = s;
                    dq_cut[ci] = gpu_dQ_[sep].data()[ci * Dim + ExpDim];
                }
            }

            cmd = [mtl_queue_ commandBuffer];
            {
                auto* descA = detail::makeDesc(B,   ExpDim);
                auto* descB = detail::makeDesc(N_s, ExpDim);
                auto* descC = detail::makeDesc(B,   N_s);

                auto* matA = detail::asMPS(gpu_dQ_se_);
                auto* matB = detail::asMPS(gpu_QS_se_);
                auto* matC = detail::asMPS(gpu_fm_);

                auto* gemm = [[MPSMatrixMultiplication alloc]
                    initWithDevice:     mtl_device_
                    transposeLeft:      NO
                    transposeRight:     YES
                    resultRows:         B
                    resultColumns:      N_s
                    interiorColumns:    ExpDim
                    alpha:              1.0
                    beta:               0.0];

                [gemm encodeToCommandBuffer: cmd
                               leftMatrix:  matA
                              rightMatrix:  matB
                             resultMatrix:  matC];

                id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
                detail::encodeElemwise2D(enc, ps_col_sub_,
                                 gpu_fm_.buf, gpu_self_.buf,
                                 static_cast<uint32_t>(B),
                                 static_cast<uint32_t>(N_s));

                [enc setComputePipelineState: ps_fm_fma_];
                [enc setBuffer: gpu_fm_.buf    offset: 0 atIndex: 0];
                [enc setBuffer: gpu_K_.buf     offset: 0 atIndex: 1];
                [enc setBuffer: gpu_K_ncs_.buf offset: 0 atIndex: 2];
                [enc setBuffer: gpu_dQ_cut_.buf offset: 0 atIndex: 3];
                uint32_t ns32 = static_cast<uint32_t>(N_s);
                [enc setBytes: &ns32 length: sizeof(uint32_t) atIndex: 4];
                const MTLSize grid = MTLSizeMake(N_s, B, 1);
                const MTLSize tpg  = MTLSizeMake(
                    std::min((NSUInteger)16, (NSUInteger)N_s),
                    std::min((NSUInteger)16, (NSUInteger)B), 1);
                [enc dispatchThreads: grid threadsPerThreadgroup: tpg];
                [enc endEncoding];
            }
            [cmd commit];
            [cmd waitUntilCompleted];

            const float* fm = gpu_fm_.data();
            for (size_t ci = 0; ci < B; ++ci) {
                const size_t atom_i  = clusters[batch_start + ci].atom_indexes[ai];
                const size_t atom_j  = clusters[batch_start + ci].atom_indexes[aj];
                const auto&  sep_geo = clusters[batch_start + ci].separations[sep];
                for (size_t si = 0; si < N_s; ++si) {
                    const Real fmv = static_cast<Real>(fm[ci * N_s + si]);
                    result.force(si, atom_i) += sep_geo.direction * fmv;
                    result.force(si, atom_j) -= sep_geo.direction * fmv;
                    result.virials(si)        += sep_geo.virials   * fmv;
                }
            }
        }
    }
}

template class BatchedNBodyGapComponent<4, 3, HasCentralAtom, 3, 1>;

} // namespace jgap
#endif // JGAP_USE_MPS