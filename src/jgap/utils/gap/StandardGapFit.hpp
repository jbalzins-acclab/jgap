#ifndef JGAP_STANDARDGAPFIT_HPP
#define JGAP_STANDARDGAPFIT_HPP

#include "StandardGapParams.hpp"
#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/core/cutoff/CosCutoff.hpp"
#include "jgap/core/kernels/SquaredExpKernel.hpp"
#include "jgap/core/potentials/CompositePotential.hpp"
#include "jgap/core/potentials/coulomb/ScreenedCoulombPotential.hpp"
#include "jgap/core/potentials/gap/GapPotential.hpp"
#include "jgap/core/potentials/gap/component/ManyBodyGapComponent.hpp"
#include "jgap/core/potentials/gap/component/ThreeBodyGapComponent.hpp"
#include "jgap/core/potentials/gap/component/TwoBodyGapComponent.hpp"
#include "jgap/core/potentials/isolated/IsolatedAtomPotential.hpp"
#include "jgap/core/sparsification/HistogramUniformSparsifier.hpp"
#include "jgap/core/transform/nbody/2b/PairDistanceTransformation.hpp"
#include "jgap/core/transform/nbody/2b/eam/CoscutoffPairFunction.hpp"
#include "jgap/core/transform/nbody/2b/eam/FSGenPairFunction.hpp"
#include "jgap/core/transform/nbody/2b/eam/PolycutoffPairFunction.hpp"
#include "jgap/core/transform/nbody/3b/Angle3bTransformation.hpp"
#include "jgap/experimental/fit/gap/StreamingQrGapFit.hpp"
#include "jgap/ext/fit/gap/QRGapFit.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"
#include "jgap/utils/gap/GapComponentUtils.hpp"

namespace jgap::utils {

    inline ValuePtr<EamPairFunction> makeStandardEamPairFunction(EamPairFunctionType type, Real cutoff2) {
        switch (type) {
            case EamPairFunctionType::FSGen2:
                return FSGenPairFunction(cutoff2, 2.0_r);
            case EamPairFunctionType::FSGen3:
                return FSGenPairFunction(cutoff2, 3.0_r);
            case EamPairFunctionType::Coscutoff:
                return CoscutoffPairFunction(cutoff2, 0.0_r);
            case EamPairFunctionType::Polycutoff:
                return PolycutoffPairFunction(cutoff2, 0.0_r);
        }
        JGAP_LOG_AND_THROW("Unknown EamPairFunctionType");
    }

    inline void standardGapFit(
        const std::string& filename,
        const std::vector<Atoms>& training_data,
        const std::vector<Regularization>& sigmas,
        const StandardGapParams& params = {}
    ) {
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
                createTwoBodyComponents<2, SquaredExpKernel<1, 1>>(training_data, trans2, kernel2, sparsifier2)
            );
        }

        // ====================================================================================
        // ManyBodyGapComponent with EAM Pair Function
        // ====================================================================================
        if (params.eam_n_sparse > 0) {
            auto eam_pf = makeStandardEamPairFunction(params.eam_pair_function, params.cutoff2);
            auto kernel_eam = SquaredExpKernel<1, 0>(1.0_r, {1.0_r});
            auto sparsifier_eam = HistogramUniformSparsifier<1>(
                params.seed, params.eam_n_sparse, std::nullopt, std::nullopt, Descriptor<1>{params.eam_min_density}
            );
            potential.addComponents(
                createEamComponents<SquaredExpKernel<1, 0>>(
                    eam_pf, kernel_eam, sparsifier_eam, training_data, params.eam_mode
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
                createThreeBodyComponents<4, SquaredExpKernel<3, 1>>(training_data, trans3, kernel3, sparsifier3)
            );
        }

        IsolatedAtomPotential isolated_atom_pot{training_data};
        ScreenedCoulombPotential sc_pot =
            params.screened_coulomb_dataset_file
                ? ScreenedCoulombPotential{*params.screened_coulomb_dataset_file, training_data}
                : ScreenedCoulombPotential{training_data};

        CompositePotential external{{
            {"isolated", isolated_atom_pot},
            {"screened_coulomb", sc_pot},
        }};

        potential.optional_external_potential = external;

        // RAM estimation & solver selection
        size_t M = 0;
        for (const auto& comp: potential.getComponents()) {
            M += comp->nSparsePoints();
        }

        size_t N = 0;
        for (const auto& atoms: training_data) {
            if (atoms.getEnergy().has_value()) {
                N += 1;
            }
            if (atoms.getForces().has_value()) {
                N += 3 * atoms.nAtoms();
            }
            if (atoms.getVirials().has_value()) {
                N += 6;
            }
        }

        if (params.approx_ram_limit_gb.has_value()) {
            const double max_bytes = *params.approx_ram_limit_gb * 1024.0 * 1024.0 * 1024.0;
            const double mm_bytes = static_cast<double>(M) * static_cast<double>(M) * sizeof(Real);
            const double two_mm_bytes = 2.0 * mm_bytes;
            const double nm_bytes = static_cast<double>(N + M) * static_cast<double>(M) * sizeof(Real);

            if (mm_bytes > max_bytes) {
                JGAP_LOG_AND_THROW("Cannot fit within the RAM requirements");
            }
            if (nm_bytes <= max_bytes) {
                QRGapFit fitter;
                fitter.fit(potential, training_data, sigmas);
            } else {
                const double remaining_rows =
                    (max_bytes / (static_cast<double>(M) * sizeof(Real))) - static_cast<double>(M);
                const size_t target_chunk_rows =
                    std::max<size_t>(1, std::min<size_t>(N, static_cast<size_t>(remaining_rows / 2.0)));

                JGAP_LOG_INFO(
                    "Using StreamingQrGapFit with target chunk rows = {} for approx_ram_limit_gb = {:.4f} GB",
                    target_chunk_rows,
                    *params.approx_ram_limit_gb
                );
                StreamingQrGapFit fitter(1e-8, target_chunk_rows);
                fitter.fit(potential, training_data, sigmas);
            }
        } else {
            QRGapFit fitter;
            fitter.fit(potential, training_data, sigmas);
        }

        SerializationRegistry<Potential>::serialize(potential, filename);
        JGAP_LOG_INFO("Saved fitted potential to {}", filename);
    }
}

#endif
