//
// Created by Jegors Balzins on 15.4.2026.
//

#include "Angle3bTransformer.hpp"

namespace jgap {

    using ToBeCalculatedPerDimension = SeparationTransformer::FixedDependencies<3, 4>::ToBeCalculatedPerDimension;

    std::array<ToBeCalculatedPerDimension, 3 + 1> Angle3bTransformer::evaluate(const TransformFrom &separations) {

        double r01 = separations.between(0, 1).magnitude;
        double r02 = separations.between(0, 2).magnitude;
        double r12 = separations.between(1, 2).magnitude;

        double f_cut_01 = cutoff->evaluate(r01), df_cut_01 = cutoff->differentiate(r01);
        double f_cut_02 = cutoff->evaluate(r02), df_cut_02 = cutoff->differentiate(r02);

        return std::array{
            ToBeCalculatedPerDimension{
                .value = r01 + r02,
                .derivatives = {
                    1.0, // wrt r_01
                    1.0, // wrt r_02
                    0.0 // wrt r_12
                }
            },
            ToBeCalculatedPerDimension{
                .value = (r01 - r02) * (r01 - r02),
                .derivatives = {
                    2.0 * (r01 - r02), // wrt r_01
                    1.0, // wrt r_02
                    0.0 // wrt r_12
                },
            },
            ToBeCalculatedPerDimension{
                .value = r12,
                .derivatives = {
                    0.0, // wrt r_01
                    0.0, // wrt r_02
                    1.0 // wrt r_12
                }
            },
            ToBeCalculatedPerDimension{
                .value = f_cut_01 * f_cut_02,
                .derivatives = {
                    df_cut_01 * f_cut_02, // wrt r_01
                    df_cut_02 * f_cut_01, // wrt r_02
                    0.0 // wrt r_12
                }
            }
        };
    }
}
