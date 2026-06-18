#ifndef JGAP_STANDARDGAPFIT_HPP
#define JGAP_STANDARDGAPFIT_HPP

#include "StandardGapParams.hpp"
#include "core/atomic/Atoms.hpp"
#include "core/cutoff/CosCutoff.hpp"
#include "core/potentials/gap/GapPotential.hpp"
#include "core/potentials/gap/component/ManyBodyGapComponent.hpp"
#include "core/potentials/gap/component/NBodyGapComponent.hpp"
#include "core/transform/eam/FSGenPairFunction.hpp"
#include "core/transform/3b/Angle3bTransformation.hpp"
#include "core/transform/2b/TwoBodyTransformation.hpp"
#include "core/kernels/SquaredExpKernel.hpp"
#include "core/sparsification/HistogramUniformSparsifier.hpp"
#include "core/fit/gap/QRGapFit.hpp"
#include "core/fit/gap/regularization/SimpleRegularizationRules.hpp"
#include "core/potentials/CompositePotential.hpp"
#include "core/potentials/zbl/ZblPotential.hpp"

namespace jgap {
    inline GapPotential standardGapFit(const std::vector<Atoms>& training_data, const StandardGapParams& params) {
        // ====================================================================================
        // ManyBodyGapComponent with EAM Pair Function
        // ====================================================================================
        auto kernel_eam = SquaredExpKernel<1, 0>(1.0, {1.0});
        auto sparsifier_eam = HistogramUniformSparsifier<1>(params.seed, params.eam_n_sparse);
        auto eam_components = EamPairFunction::createComponents(
            params.eam_pf, kernel_eam, sparsifier_eam, training_data, params.eam_mode);

        // ====================================================================================
        // 3-Body and 2-Body Components
        // ====================================================================================
        auto trans3 = Angle3bTransformation(CosCutoff(params.cutoff3, params.cutoff3_width));
        auto kernel3 = SquaredExpKernel<3, 1>(1.0, {1.0, 1.0, 1.0});
        auto sparsifier3 = HistogramUniformSparsifier<4>(params.seed, params.n_sparse3, params.active3);
        auto components3 = NBodyGapComponent<4, 3, HasCentralAtom, SquaredExpKernel<3, 1>>::createComponents(
            training_data, trans3, kernel3, sparsifier3);

        auto trans2 = TwoBodyTransformation(CosCutoff(params.cutoff2, params.cutoff2_width));
        auto kernel2 = SquaredExpKernel<1, 1>(10.0, {1.0});
        auto sparsifier2 = HistogramUniformSparsifier<2>(params.seed, params.n_sparse2, params.active2);
        auto components2 = NBodyGapComponent<2, 2, Symmetric, SquaredExpKernel<1, 1>>::createComponents(
            training_data, trans2, kernel2, sparsifier2);

        // ====================================================================================
        // Construct Potential
        // ====================================================================================
        GapPotential potential;
        potential.addComponents(eam_components);
        potential.addComponents(components3);
        potential.addComponents(components2);

        IsolatedAtomPotential isolated_atom_pot{training_data, false};
        ZblPotential zbp_pot{training_data};

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