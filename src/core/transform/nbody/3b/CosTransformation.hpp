#ifndef JGAP_COSTRANSFORMATION_HPP
#define JGAP_COSTRANSFORMATION_HPP
#include "ThreeBodyTransformation.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {
    class CosTransformation final : public ThreeBodyTransformation<3> {
    public:
        Descriptor<3> evaluate(const Cluster3& cluster) const override {
            Real r01 = cluster.separationBetween(0, 1);
            Real r02 = cluster.separationBetween(0, 2);
            Real r12 = cluster.separationBetween(1, 2);

            Real cos12 = (r01 * r01 + r02 * r02 - r12 * r12) / (2.0 * r01 * r02);

            return {r01, r02, cos12};
        }

        ThreeBodyDescriptor<3> evaluateAndDifferentiate(const Cluster3& cluster) const override {
            Real r01 = cluster.separationBetween(0, 1);
            Real r02 = cluster.separationBetween(0, 2);
            Real r12 = cluster.separationBetween(1, 2);

            Real inv_r01 = 1.0 / r01;
            Real inv_r02 = 1.0 / r02;
            Real cos12 = (r01 * r01 + r02 * r02 - r12 * r12) * (0.5 * inv_r01 * inv_r02);

            // clang-format off
            return {
                .value = {
                    r01,
                    r02,
                    cos12
                },
                .derivatives = {
                    std::array{ // wrt r_01
                        1.0,
                        0.0,
                        inv_r02 - cos12 * inv_r01
                    },
                    std::array{ // wrt r_02
                        0.0,
                        1.0,
                        inv_r01 - cos12 * inv_r02
                    },
                    std::array{ // wrt r_12
                        0.0,
                        0.0,
                        -r12 * inv_r01 * inv_r02
                    }
                }
            };
            // clang-format on
        }

        Cutoffs getCutoffs() const override {
            JGAP_LOG_AND_THROW("Cutoff is implicit; note: CosTransformation is intended for TabGap components only");
        }

        CosTransformation* clone() const override { return new CosTransformation(*this); }

        bool isSwapInvariant(size_t idx1, size_t idx2) const override { return false; }
    };
}
#endif
