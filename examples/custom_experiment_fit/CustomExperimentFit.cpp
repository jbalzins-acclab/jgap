// Example: define a brand new custom kernel (FractionalExpKernel with power 1.5) in user code,
// fit a 2-body only GAP using ElementIncrementalQRGapFit, do not save the GAP potential,
// but tabulate directly and save the tabulated potential (.tabgap.h5).
//
// Usage: custom_experiment_fit <training.xyz> <output_prefix> [--ram-limit <gb>]
//   writes <output_prefix>.tabgap.h5

#include <chrono>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "jgap/core/atomic/species/composition/Species2Sorted.hpp"
#include "jgap/core/fit/gap/regularization/PerConfigTypeRegularizationRules.hpp"
#include "jgap/core/kernels/Kernel.hpp"
#include "jgap/core/potentials/gap/component/TwoBodyGapComponent.hpp"
#include "jgap/impl/cutoff/CosCutoff.hpp"
#include "jgap/impl/fit/gap/ElementIncrementalQRGapFit.hpp"
#include "jgap/impl/transform/nbody/2b/PairDistanceTransformation.hpp"
#include "jgap/jgap.hpp"

using namespace jgap;
using namespace jgap::utils;

// =============================================================================
// User-Defined Custom Kernel: Fractional Exponential Kernel (power p = 1.5)
//   k(x, y) = \sigma^2 * exp( - ||(x - y) / \ell||^1.5 ) * cutoff(x) * cutoff(y)
//
// By Schoenberg's theorem (1938), exp(-||r||^p) is strictly positive definite
// if and only if 0 < p <= 2. The power p = 1.5 guarantees positive definiteness.
// =============================================================================
template<size_t ExpDimensions, size_t CutoffDimensions>
    requires(CutoffDimensions <= 1)
class FractionalExpKernel final : public Kernel<ExpDimensions + CutoffDimensions> {
public:
    static constexpr size_t ExpDim = ExpDimensions;
    static constexpr size_t CutoffDim = CutoffDimensions;
    static constexpr size_t TotalDimensions = ExpDimensions + CutoffDimensions;

    using KernelValueAndGradient = typename Kernel<TotalDimensions>::KernelValueAndGradient;

    FractionalExpKernel() = default;

    FractionalExpKernel(const double energy_scale, const std::array<double, ExpDimensions>& length_scales) {
        prefactor = energy_scale * energy_scale;
        for (size_t dim = 0; dim < ExpDimensions; dim++) {
            inverse_length_scales_squared[dim] = 1.0 / (length_scales[dim] * length_scales[dim]);
        }
    }

    double getEnergyScale() const { return std::sqrt(prefactor); }

    double value(const Descriptor<TotalDimensions>& q1, const Descriptor<TotalDimensions>& q2) const override {
        return Kernel<TotalDimensions>::value(q1, q2);
    }

    KernelValueAndGradient valueAndGradient(
        const Descriptor<TotalDimensions>& sparse_point, const Descriptor<TotalDimensions>& q
    ) const override {
        double dist_sq = 0.0;
        for (size_t dim = 0; dim < ExpDimensions; dim++) {
            double diff = q[dim] - sparse_point[dim];
            dist_sq += diff * diff * inverse_length_scales_squared[dim];
        }
        double u = std::sqrt(dist_sq);
        double base_val = prefactor * std::exp(-std::pow(u, 1.5));
        double val = base_val;

        std::array<double, TotalDimensions> gradient{};
        if constexpr (CutoffDimensions == 1) {
            gradient[ExpDimensions] = base_val * sparse_point[ExpDimensions];
            val = base_val * sparse_point[ExpDimensions] * q[ExpDimensions];
        }

        if (u > 1e-12) {
            double factor = 1.5 * val / std::sqrt(u);
            for (size_t dim = 0; dim < ExpDimensions; dim++) {
                gradient[dim] = factor * (sparse_point[dim] - q[dim]) * inverse_length_scales_squared[dim];
            }
        }

        return {.value = val, .gradient = gradient};
    }

    FractionalExpKernel* clone() const override { return new FractionalExpKernel(*this); }

private:
    double prefactor{1.0};
    std::array<double, ExpDimensions> inverse_length_scales_squared{};
};

int main(int argc, char** argv) {
    CurrentLogger::initDefault({.stdout_log_debug = true});

    std::vector<std::string> positional_args;
    double ram_limit_gb = 2.0;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--ram-limit" && i + 1 < argc) {
            ram_limit_gb = std::stod(argv[++i]);
        } else if (arg.starts_with("--ram-limit=")) {
            ram_limit_gb = std::stod(arg.substr(12));
        } else if (arg.starts_with("ram_limit=")) {
            ram_limit_gb = std::stod(arg.substr(10));
        } else if (arg == "-h" || arg == "--help") {
            std::cout << "Usage: " << argv[0] << " <training.xyz> <output_prefix> [--ram-limit <gb>]\n"
                      << "  --ram-limit <gb>   RAM limit in GB for ElementIncrementalQRGapFit (default: 2.0)\n";
            return 0;
        } else {
            positional_args.push_back(arg);
        }
    }

    if (positional_args.size() != 2) {
        std::cerr << "Usage: " << argv[0] << " <training.xyz> <output_prefix> [--ram-limit <gb>]\n";
        return 1;
    }
    const std::string training_file = positional_args[0];
    const std::string output_prefix = positional_args[1];

    const auto start = std::chrono::steady_clock::now();

    JGAP_LOG_INFO("Custom Experiment 2b-only fit with FractionalExpKernel (-|r/l|^1.5) on {} (RAM limit: {} GB)",
                  training_file, ram_limit_gb);
    auto training_data = Atoms::readAtoms(training_file);

    const size_t seed = 42;
    const double cutoff2 = 5.0;
    const double cutoff2_width = 1.0;
    const size_t n_sparse2 = 15;

    GapPotential potential;

    // 2-Body Components with Custom FractionalExpKernel
    auto trans2 = PairDistanceTransformation(CosCutoff(cutoff2, cutoff2_width));
    auto custom_kernel = FractionalExpKernel<1, 1>(1.0, {1.5});
    auto sparsifier2 = HistogramUniformSparsifier<2>(seed, n_sparse2, std::array{true, false});

    std::set<Species2Sorted> pairs;
    for (const auto& atoms : training_data) {
        NeighbourLists nl(atoms, cutoff2);
        auto sets = Species2Sorted::getAll(nl);
        pairs.insert(sets.begin(), sets.end());
    }
    for (const auto& pair : pairs) {
        potential.addComponent(TwoBodyGapComponent<2, FractionalExpKernel<1, 1>>(
            pair, trans2, custom_kernel, sparsifier2, training_data
        ));
    }

    // Isolated atom external potential
    potential.optional_external_potential = IsolatedAtomPotential{training_data};

    // Fit with ElementIncrementalQRGapFit
    ElementIncrementalQRGapFit fitter(ram_limit_gb);
    PerConfigTypeRegularizationRules rules{PerConfigTypeSigmas(0.001, 0.05, 0.1, 0.02)};
    auto sigmas = rules.determineForAll(training_data);
    fitter.fit(potential, training_data, sigmas);

    // Note: We do NOT serialize the GAP potential here, only tabulate and save the tabGAP potential
    JGAP_LOG_INFO("Skipping GAP serialization; tabulating directly to .tabgap.h5...");
    StandardTabulationParams tab_params;
    tab_params.n_grid_2b = 5000;
    standardTabulation(potential, output_prefix, tab_params);

    std::cout << "Execution time: " << formatDuration(elapsedMillisSince(start)) << std::endl;
    return 0;
}
