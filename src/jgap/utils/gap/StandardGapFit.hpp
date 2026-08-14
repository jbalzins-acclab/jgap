#ifndef JGAP_STANDARDGAPFIT_HPP
#define JGAP_STANDARDGAPFIT_HPP

#include "StandardGapParams.hpp"
#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/core/cutoff/CosCutoff.hpp"
#include "jgap/core/kernels/SquaredExpKernel.hpp"
#include "jgap/core/potentials/CompositePotential.hpp"
#include "jgap/core/potentials/gap/GapPotential.hpp"
#include "jgap/core/potentials/gap/component/AtomicThreeBodyGapComponent.hpp"
#include "jgap/core/potentials/gap/component/ManyBodyGapComponent.hpp"
#include "jgap/core/potentials/gap/component/TwoBodyGapComponent.hpp"
#include "jgap/core/potentials/isolated/IsolatedAtomPotential.hpp"
#include "jgap/core/potentials/coulomb/ScreenedCoulombPotential.hpp"
#include "jgap/core/sparsification/HistogramUniformSparsifier.hpp"
#include "jgap/core/transform/nbody/2b/PairDistanceTransformation.hpp"
#include "jgap/core/transform/nbody/3b/Angle3bTransformation.hpp"
#include "jgap/ext/fit/gap/QRGapFit.hpp"
#include "jgap/utils/gap/GapComponentUtils.hpp"

namespace jgap {
    inline GapPotential standardGapFit(const std::vector<Atoms>& training_data, const StandardGapParams& params) {
        GapPotential potential;

        if (params.n_sparse2 == 0 && params.n_sparse3 == 0 && params.eam_n_sparse == 0) {
            JGAP_LOG_AND_THROW("Cannot make a standard GAP potential without any components");
        }

        // ====================================================================================
        // 2-Body Components
        // ====================================================================================
        if (params.n_sparse2 > 0) {
            auto trans2 = PairDistanceTransformation(CosCutoff(params.cutoff2, params.cutoff2_width));
            auto kernel2 = SquaredExpKernel<1, 1>(10.0_r, {1.0_r});
            auto sparsifier2 = HistogramUniformSparsifier<2>(params.seed, params.n_sparse2, std::array{true, false});
            potential.addComponents(
                createTwoBodyComponents<2, SquaredExpKernel<1, 1>>(
                    training_data, trans2, kernel2, sparsifier2
                )
            );
        }

        // ====================================================================================
        // ManyBodyGapComponent with EAM Pair Function
        // ====================================================================================
        if (params.eam_n_sparse > 0) {
            auto kernel_eam = SquaredExpKernel<1, 0>(1.0_r, {1.0_r});
            auto sparsifier_eam = HistogramUniformSparsifier<1>(params.seed, params.eam_n_sparse);
            potential.addComponents(
                createEamComponents<SquaredExpKernel<1, 0>>(
                    params.eam_pf, kernel_eam, sparsifier_eam, training_data, params.eam_mode
                )
            );
        }

        // ====================================================================================
        // 3-Body Components
        // ====================================================================================
        if (params.n_sparse3 > 0) {
            auto trans3 = Angle3bTransformation(CosCutoff(params.cutoff3, params.cutoff3_width));
            auto kernel3 = SquaredExpKernel<3, 1>(1.0_r, {1.0_r, 1.0_r, 1.0_r});
            auto sparsifier3 =
                HistogramUniformSparsifier<4>(params.seed, params.n_sparse3, std::array{true, true, true, false});
            potential.addComponents(
                createThreeBodyComponents<4, SquaredExpKernel<3, 1>>(
                    training_data, trans3, kernel3, sparsifier3
                )
            );
        }

        IsolatedAtomPotential isolated_atom_pot{training_data};
        ScreenedCoulombPotential sc_pot = params.screened_coulomb_dataset_file
                                             ? ScreenedCoulombPotential{*params.screened_coulomb_dataset_file, training_data}
                                             : ScreenedCoulombPotential{training_data};

        CompositePotential external{{
            {"isolated", isolated_atom_pot},
            {"screened_coulomb", sc_pot},
        }};

        potential.optional_external_potential = external;

        // Fit
        QRGapFit fitter;
        fitter.fit(potential, training_data, *params.regularization_rules);

        return potential;
    }
}

#endif
