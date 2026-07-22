#ifndef JGAP_STANDARDGAPFIT_HPP
#define JGAP_STANDARDGAPFIT_HPP

#include "StandardGapParams.hpp"
#include "core/atomic/Atoms.hpp"
#include "core/cutoff/CosCutoff.hpp"
#include "core/kernels/SquaredExpKernel.hpp"
#include "core/potentials/CompositePotential.hpp"
#include "core/potentials/gap/GapPotential.hpp"
#include "core/potentials/gap/component/AtomicThreeBodyGapComponent.hpp"
#include "core/potentials/gap/component/ManyBodyGapComponent.hpp"
#include "core/potentials/gap/component/TwoBodyGapComponent.hpp"
#include "core/potentials/isolated/IsolatedAtomPotential.hpp"
#include "core/potentials/zbl/ZblPotential.hpp"
#include "core/sparsification/HistogramUniformSparsifier.hpp"
#include "core/transform/nbody/2b/PairDistanceTransformation.hpp"
#include "core/transform/nbody/3b/Angle3bTransformation.hpp"
#include "ext/fit/gap/QRGapFit.hpp"

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
            auto kernel2 = SquaredExpKernel<1, 1>(10.0, {1.0});
            auto sparsifier2 = HistogramUniformSparsifier<2>(params.seed, params.n_sparse2, std::array{true, false});
            potential.addComponents(TwoBodyGapComponent<2, SquaredExpKernel<1, 1>>::createComponents(
                training_data, trans2, kernel2, sparsifier2));
        }

        // ====================================================================================
        // ManyBodyGapComponent with EAM Pair Function
        // ====================================================================================
        if (params.eam_n_sparse > 0) {
            auto kernel_eam = SquaredExpKernel<1, 0>(1.0, {1.0});
            auto sparsifier_eam = HistogramUniformSparsifier<1>(params.seed, params.eam_n_sparse);
            potential.addComponents(EamPairFunction::createComponents(params.eam_pf, kernel_eam, sparsifier_eam,
                                                                      training_data, params.eam_mode));
        }

        // ====================================================================================
        // 3-Body Components
        // ====================================================================================
        if (params.n_sparse3 > 0) {
            auto trans3 = Angle3bTransformation(CosCutoff(params.cutoff3, params.cutoff3_width));
            auto kernel3 = SquaredExpKernel<3, 1>(1.0, {1.0, 1.0, 1.0});
            auto sparsifier3 =
                HistogramUniformSparsifier<4>(params.seed, params.n_sparse3, std::array{true, true, true, false});
            potential.addComponents(AtomicThreeBodyGapComponent<4, SquaredExpKernel<3, 1>>::createComponents(
                training_data, trans3, kernel3, sparsifier3));
        }

        IsolatedAtomPotential isolated_atom_pot{training_data};
        ZblPotential zbp_pot = params.zbl_dataset_file ? ZblPotential{*params.zbl_dataset_file, training_data}
                                                       : ZblPotential{training_data};

        CompositePotential external{{
            {"isolated", isolated_atom_pot},
            {"zbl", zbp_pot},
        }};

        potential.optional_external_potential = external;

        // Fit
        QRGapFit fitter;
        fitter.fit(potential, training_data, *params.regularization_rules);

        return potential;
    }
}

#endif
