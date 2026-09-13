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
            double r01 = cluster.separation01.magnitude;
            double r02 = cluster.separation02.magnitude;
            double r12 = cluster.separation12.magnitude;

            double inv_r01 = 1.0 / r01;
            double inv_r02 = 1.0 / r02;
            double cos12 = (r01 * r01 + r02 * r02 - r12 * r12) * (0.5 * inv_r01 * inv_r02);

            const auto& dir01 = cluster.separation01.direction;
            const auto& dir02 = cluster.separation02.direction;
            const auto& dir12 = cluster.separation12.direction;

            double dcos_dr01 = inv_r02 - cos12 * inv_r01;
            double dcos_dr02 = inv_r01 - cos12 * inv_r02;
            double dcos_dr12 = -r12 * inv_r01 * inv_r02;

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
