#ifndef JGAP_COSTRANSFORMATION_HPP
#define JGAP_COSTRANSFORMATION_HPP
#include "../../../io/log/CurrentLogger.hpp"
#include "ThreeBodyTransformation.hpp"

namespace jgap {
    class CosTransformation final : public ThreeBodyTransformation<3> {
    public:
        Descriptor<3> evaluate(const Cluster3& cluster) const override {
            return ThreeBodyTransformation<3>::evaluate(cluster);
        }

        ThreeBodyDescriptor<3> evaluateAndDifferentiate(const Cluster3& cluster) const override final {
            Real r01 = cluster.separation01.magnitude;
            Real r02 = cluster.separation02.magnitude;
            Real r12 = cluster.separation12.magnitude;

            Real inv_r01 = 1.0_r / r01;
            Real inv_r02 = 1.0_r / r02;
            Real cos12 = (r01 * r01 + r02 * r02 - r12 * r12) * (0.5_r * inv_r01 * inv_r02);

            const auto& dir01 = cluster.separation01.direction;
            const auto& dir02 = cluster.separation02.direction;
            const auto& dir12 = cluster.separation12.direction;

            Real dcos_dr01 = inv_r02 - cos12 * inv_r01;
            Real dcos_dr02 = inv_r01 - cos12 * inv_r02;
            Real dcos_dr12 = -r12 * inv_r01 * inv_r02;

            return { 
                .value = { 
                    r01, 
                    r02, 
                    cos12,
                },
                .grad_r1 = { 
                    dir01, 
                    Vector3{}, 
                    dcos_dr01 * dir01 - dcos_dr12 * dir12
                },
                .grad_r2 = { 
                    Vector3{}, 
                    dir02, 
                    dcos_dr02 * dir02 + dcos_dr12 * dir12
                } 
            };
        }

        Cutoffs getCutoffs() const override {
            JGAP_LOG_AND_THROW("Cutoff is implicit; note: CosTransformation is intended for TabGap components only");
        }
        bool isRotationallyInvariant() const override { return true; }

        CosTransformation* clone() const override { return new CosTransformation(*this); }

        bool isSwapInvariant(size_t idx1, size_t idx2) const override { return false; }
    };
}
#endif
