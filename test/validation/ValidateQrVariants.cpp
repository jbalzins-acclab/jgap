#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "common/ValidationConfig.hpp"
#include "common/ValidationUtils.hpp"

#include "jgap/core/UnseqFor.hpp"
#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/core/atomic/iteration/Cluster2Expansion.hpp"
#include "jgap/core/atomic/iteration/Cluster3Expansion.hpp"
#include "jgap/core/atomic/neighbours/NeighbourLists.hpp"
#include "jgap/impl/cutoff/CosCutoff.hpp"
#include "jgap/core/fit/gap/regularization/PerConfigTypeRegularizationRules.hpp"
#include "jgap/impl/kernels/SquaredExpKernel.hpp"
#include "jgap/core/potentials/CompositePotential.hpp"
#include "jgap/core/potentials/gap/GapPotential.hpp"
#include "jgap/core/potentials/gap/component/ManyBodyGapComponent.hpp"
#include "jgap/core/potentials/gap/component/ThreeBodyGapComponent.hpp"
#include "jgap/core/potentials/gap/component/TwoBodyGapComponent.hpp"
#include "jgap/core/potentials/isolated/IsolatedAtomPotential.hpp"
#include "jgap/core/sparsification/HistogramUniformSparsifier.hpp"
#include "jgap/impl/transform/nbody/2b/PairDistanceTransformation.hpp"
#include "jgap/impl/transform/nbody/2b/eam/PolycutoffPairFunction.hpp"
#include "jgap/impl/transform/nbody/3b/Angle3bTransformation.hpp"
#include "jgap/impl/fit/gap/BlockIncrementalQRGapFit.hpp"
#include "jgap/impl/fit/gap/ElementIncrementalQRGapFit.hpp"
#include "jgap/impl/fit/gap/QRGapFit.hpp"
#include "jgap/io/convert/QuipXmlConverter.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"
#include "jgap/utils/gap/GapComponentUtils.hpp"

using namespace jgap;
using namespace jgap::validation;
namespace fs = std::filesystem;

struct PrecomputedDescriptors {
    std::vector<Species2Sorted> pairs;
    std::vector<Species> eam_species;
    std::vector<Species3AtomicSorted> triplets;

    std::vector<std::vector<std::vector<Descriptor<2>>>> d2_per_frame;
    std::vector<std::vector<std::vector<Descriptor<1>>>> deam_per_frame;
    std::vector<std::vector<std::vector<Descriptor<4>>>> d3_per_frame;

    PairDistanceTransformation trans2{CosCutoff(5.0, 1.0)};
    SquaredExpKernel<1, 1> kernel2{10.0, {1.0}};

    PolycutoffPairFunction eam_pf{5.0, 0.0};
    SquaredExpKernel<1, 0> kernel_eam{1.0, {1.0}};

    Angle3bTransformation trans3{CosCutoff(4.0, 0.6)};
    SquaredExpKernel<3, 1> kernel3{1.0, {1.0, 1.0, 1.0}};

    std::map<Species, ValuePtr<NBodyAggregator<1>>> eam_aggs;
};

static PrecomputedDescriptors precomputeAllDescriptors(
    const std::vector<Atoms>& all_atoms, const std::vector<std::optional<NeighbourLists>>& nls
) {
    PrecomputedDescriptors p;
    p.pairs = {Species2Sorted("Fe", "Fe"), Species2Sorted("Fe", "Ni"), Species2Sorted("Ni", "Ni")};
    p.eam_species = {Species("Fe"), Species("Ni")};
    p.triplets = {
        Species3AtomicSorted(Species("Fe"), Species("Fe"), Species("Fe")),
        Species3AtomicSorted(Species("Fe"), Species("Fe"), Species("Ni")),
        Species3AtomicSorted(Species("Fe"), Species("Ni"), Species("Ni")),
        Species3AtomicSorted(Species("Ni"), Species("Fe"), Species("Fe")),
        Species3AtomicSorted(Species("Ni"), Species("Fe"), Species("Ni")),
        Species3AtomicSorted(Species("Ni"), Species("Ni"), Species("Ni"))
    };

    p.eam_aggs = utils::createEamAggregators(p.eam_pf, all_atoms, EamMode::Blind);

    const size_t N = all_atoms.size();
    p.d2_per_frame.resize(N, std::vector<std::vector<Descriptor<2>>>(p.pairs.size()));
    p.deam_per_frame.resize(N, std::vector<std::vector<Descriptor<1>>>(p.eam_species.size()));
    p.d3_per_frame.resize(N, std::vector<std::vector<Descriptor<4>>>(p.triplets.size()));

    unseqForIndex(0, N, [&](size_t i) {
        const auto& nl = *nls[i];
        for (size_t k = 0; k < p.pairs.size(); ++k) {
            Cluster2Expansion exp(p.pairs[k]);
            exp.forEach(nl, [&](const Cluster2& c) { p.d2_per_frame[i][k].push_back(p.trans2.evaluate(c)); });
        }
        for (size_t k = 0; k < p.eam_species.size(); ++k) {
            auto agg = p.eam_aggs[p.eam_species[k]]->aggregate(nl);
            for (const auto& v: agg.values) {
                p.deam_per_frame[i][k].push_back(Descriptor<1>{v[0]});
            }
        }
        for (size_t k = 0; k < p.triplets.size(); ++k) {
            Cluster3Expansion exp(
                p.triplets[k],
                p.trans3.isSwapInvariant(1, 2) ? ClusterPermutationMode::NoNodePermutation
                                               : ClusterPermutationMode::PermuteSameSpeciesNodes
            );
            exp.forEach(nl, [&](const Cluster3& c) { p.d3_per_frame[i][k].push_back(p.trans3.evaluate(c)); });
        }
    });

    return p;
}

static GapPotential createFeniPotential(
    const PrecomputedDescriptors& desc,
    const std::vector<size_t>& active_indices,
    const std::vector<Atoms>& active_atoms,
    size_t n_3b,
    size_t sparse_seed,
    const ValuePtr<Potential>& glue_pairpot
) {
    double ratio = static_cast<double>(n_3b) / 500.0;
    size_t n_2b = std::max<size_t>(1, std::lround(20.0 * ratio));
    size_t n_eam = std::max<size_t>(1, std::lround(20.0 * ratio));

    GapPotential potential;

    auto sparsifier2 = HistogramUniformSparsifier<2>(sparse_seed, n_2b, std::array{true, false});
    for (size_t k = 0; k < desc.pairs.size(); ++k) {
        std::vector<Descriptor<2>> comp_descs;
        for (size_t idx: active_indices) {
            const auto& frame_d = desc.d2_per_frame[idx][k];
            comp_descs.insert(comp_descs.end(), frame_d.begin(), frame_d.end());
        }
        auto sparse_pts = sparsifier2.selectSparsePoints(comp_descs);
        potential.components.push_back(
            TwoBodyGapComponent<2, SquaredExpKernel<1, 1>>(desc.pairs[k], desc.trans2, desc.kernel2, sparse_pts)
        );
    }

    auto sparsifier_eam =
        HistogramUniformSparsifier<1>(sparse_seed, n_eam, std::nullopt, std::nullopt, Descriptor<1>{0.05});
    for (size_t k = 0; k < desc.eam_species.size(); ++k) {
        std::vector<Descriptor<1>> comp_descs;
        for (size_t idx: active_indices) {
            const auto& frame_d = desc.deam_per_frame[idx][k];
            comp_descs.insert(comp_descs.end(), frame_d.begin(), frame_d.end());
        }
        auto sparse_pts = sparsifier_eam.selectSparsePoints(comp_descs);
        potential.components.push_back(
            ManyBodyGapComponent<1, SquaredExpKernel<1, 0>>(
                desc.eam_aggs.at(desc.eam_species[k]), desc.kernel_eam, sparse_pts
            )
        );
    }

    auto sparsifier3 = HistogramUniformSparsifier<4>(sparse_seed, n_3b, std::array{true, true, true, false});
    for (size_t k = 0; k < desc.triplets.size(); ++k) {
        std::vector<Descriptor<4>> comp_descs;
        for (size_t idx: active_indices) {
            const auto& frame_d = desc.d3_per_frame[idx][k];
            comp_descs.insert(comp_descs.end(), frame_d.begin(), frame_d.end());
        }
        auto sparse_pts = sparsifier3.selectSparsePoints(comp_descs);
        potential.components.push_back(
            ThreeBodyGapComponent<4, SquaredExpKernel<3, 1>>(desc.triplets[k], desc.trans3, desc.kernel3, sparse_pts)
        );
    }

    IsolatedAtomPotential isolated_pot(active_atoms);
    std::map<std::string, ValuePtr<Potential>> ext_map;
    ext_map["isolated"] = isolated_pot;
    if (glue_pairpot) {
        ext_map["glue"] = glue_pairpot;
    }
    potential.optional_external_potential = CompositePotential(ext_map);

    return potential;
}

static double computeRamLimitGb(size_t M, double factor) {
    double total_rows = factor * static_cast<double>(M);
    double ram_bytes = total_rows * static_cast<double>(M) * sizeof(double);
    return std::max(0.005, ram_bytes / (1024.0 * 1024.0 * 1024.0));
}

static std::vector<double> getAllCoefficients(const GapPotential& pot) {
    std::vector<double> coeffs;
    for (const auto& comp: pot.getComponents()) {
        const auto& c = comp->getCoefficients();
        coeffs.insert(coeffs.end(), c.begin(), c.end());
    }
    return coeffs;
}

int main(int argc, char* argv[]) {
    QrValidationOptions opts = QrValidationOptions::parse(argc, argv);
    ValidationReporter reporter("ValidateQrVariants");

    std::cout << "======================================================================\n";
    std::cout << "JGAP Validation: QR Solver Variants Benchmark\n";
    std::cout << "  3B Sparse Points : min=" << opts.min_m3b << ", max=" << opts.max_m3b << ", step=" << opts.step_m3b
              << "\n";
    std::cout << "  Seeds            : " << opts.sparse_seeds.size() << " seeds\n";
    std::cout << "  NRMSE Threshold  : <= " << opts.max_allowed_nrmse_pct << " %\n";
    std::cout << "  Cosine Threshold : >= " << opts.min_allowed_cosine_sim << "\n";
    std::cout << "======================================================================\n";

    fs::path root = findRepoRoot();
    fs::path train_xyz = root / "test" / "resources" / "structure-databases" / "feni-train.xyz";
    fs::path glue_xml = root / "test" / "resources" / "quip_potentials" / "pairpot_Fe_Ni.xml";

    if (!fs::exists(train_xyz)) {
        std::cerr << "Error: training data not found at " << train_xyz << "\n";
        return 1;
    }

    MainXYZPropertyNames prop_names;
    prop_names.virials = "virial_fit";
    const std::vector<Atoms> all_train_atoms = Atoms::readAtoms(train_xyz.string(), prop_names);

    ValuePtr<Potential> glue_pot;
    if (fs::exists(glue_xml)) {
        glue_pot = QuipXmlConverter::transform(glue_xml.string());
    }

    std::cout << "Precomputing neighbour lists and descriptors for " << all_train_atoms.size() << " frames..."
              << std::flush;
    std::vector<std::optional<NeighbourLists>> nls(all_train_atoms.size());
    unseqForIndex(0, all_train_atoms.size(), [&](size_t i) { nls[i] = NeighbourLists(all_train_atoms[i], 5.0); });
    PrecomputedDescriptors precomputed = precomputeAllDescriptors(all_train_atoms, nls);
    std::cout << " Done.\n";

    const PerConfigTypeRegularizationRules regularization(
        PerConfigTypeSigmas(0.002, 0.1, 0.2),
        "isolated_atom:0.0001:0.04:0.04:0.0:"
        "liquid:0.01:0.5:2.0:0.0:"
        "dimer:0.01:0.5:2.0:0.0:"
        "short_range:0.01:0.5:2.0:0.0:"
        "interstitial:0.005:0.2:0.5:0.0:"
        "surface:0.005:0.2:0.5:0.0:"
        "pure_bulk:0.001:0.05:0.1:0.0:"
        "defect:0.005:0.2:0.5:0.0:"
        "bulk:0.002:0.1:0.2:0.0:"
        "alloy_bulk:0.002:0.1:0.2:0.0:"
        "amorphous:0.01:0.5:2.0:0.0"
    );

    std::ofstream csv;
    if (!opts.output_csv.empty()) {
        csv.open(opts.output_csv);
        csv << "n_sparse_3b,seed,variant,m_total,fit_time_s,nrmse_pct,max_rel_pct,cosine_sim\n";
    }

    std::vector<size_t> active_indices(all_train_atoms.size());
    std::iota(active_indices.begin(), active_indices.end(), 0);

    for (size_t n_3b = opts.min_m3b; n_3b <= opts.max_m3b; n_3b += opts.step_m3b) {
        for (size_t seed: opts.sparse_seeds) {
            std::cout << "\n[RUN] N_3b=" << n_3b << ", Seed=" << seed << "...\n";

            GapPotential base_pot =
                createFeniPotential(precomputed, active_indices, all_train_atoms, n_3b, seed, glue_pot);

            size_t M = 0;
            for (const auto& comp: base_pot.components) {
                M += comp->nSparsePoints();
            }

            auto sigmas = regularization.determineForAll(all_train_atoms);

            // 1. Reference Full QR
            std::cout << "  Fitting Reference Full QR (M=" << M << ")..." << std::flush;
            GapPotential pot_ref = base_pot;
            QRGapFit full_fit(1e-8);
            full_fit.fit(pot_ref, all_train_atoms, sigmas);
            auto coeffs_ref = getAllCoefficients(pot_ref);
            std::cout << " Done.\n";

            // 2. Block Incremental QR variants
            for (double factor: {2.0, 5.0}) {
                std::string var_name = "Block_QR_M+B=" + std::to_string(static_cast<int>(factor)) + "M";
                std::cout << "  Fitting " << var_name << "..." << std::flush;
                double ram_gb = computeRamLimitGb(M, factor);
                GapPotential pot_block = base_pot;
                BlockIncrementalQRGapFit block_fit(1e-8, ram_gb);
                block_fit.fit(pot_block, all_train_atoms, sigmas);
                auto coeffs_block = getAllCoefficients(pot_block);
                std::cout << " Done.\n";

                double nrmse = computeNrmsePct(coeffs_block, coeffs_ref);
                double sigrel = computeSigRelPct(coeffs_block, coeffs_ref);
                double csim = computeCosineSim(coeffs_block, coeffs_ref);

                std::string test_id = "N3b_" + std::to_string(n_3b) + "_" + var_name;
                reporter.assertLessOrEqual(test_id, "NRMSE", nrmse, opts.max_allowed_nrmse_pct, "%");
                reporter.assertLessOrEqual(test_id, "SigRel", sigrel, opts.max_allowed_sigrel_pct, "%");
                reporter.assertGreaterOrEqual(test_id, "CosineSim", csim, opts.min_allowed_cosine_sim);

                if (csv.is_open()) {
                    csv << n_3b << "," << seed << "," << var_name << "," << M << ",0.0," << nrmse << "," << sigrel
                        << "," << csim << "\n";
                }
            }

            // 3. Element Incremental QR
            {
                std::string var_name = "Element_QR_M+B~2M";
                std::cout << "  Fitting " << var_name << "..." << std::flush;
                double ram_gb = computeRamLimitGb(M, 2.0);
                GapPotential pot_elem = base_pot;
                ElementIncrementalQRGapFit elem_fit(1e-8, ram_gb);
                elem_fit.fit(pot_elem, all_train_atoms, sigmas);
                auto coeffs_elem = getAllCoefficients(pot_elem);
                std::cout << " Done.\n";

                double nrmse = computeNrmsePct(coeffs_elem, coeffs_ref);
                double sigrel = computeSigRelPct(coeffs_elem, coeffs_ref);
                double csim = computeCosineSim(coeffs_elem, coeffs_ref);

                std::string test_id = "N3b_" + std::to_string(n_3b) + "_" + var_name;
                reporter.assertLessOrEqual(test_id, "NRMSE", nrmse, opts.max_allowed_nrmse_pct, "%");
                reporter.assertLessOrEqual(test_id, "SigRel", sigrel, opts.max_allowed_sigrel_pct, "%");
                reporter.assertGreaterOrEqual(test_id, "CosineSim", csim, opts.min_allowed_cosine_sim);

                if (csv.is_open()) {
                    csv << n_3b << "," << seed << "," << var_name << "," << M << ",0.0," << nrmse << "," << sigrel
                        << "," << csim << "\n";
                }
            }
        }
    }

    return reporter.summarizeAndExit();
}
