#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/core/fit/gap/regularization/PerConfigTypeRegularizationRules.hpp"
#include "jgap/core/potentials/gap/GapPotential.hpp"
#include "jgap/core/potentials/tabgap/TabGapPotential.hpp"
#include "jgap/experimental/fit/gap/ElementalQRGapFit.hpp"
#include "jgap/experimental/fit/gap/BlockIncrementalQRGapFit.hpp"
#include "jgap/io/convert/QuipXmlConverter.hpp"
#include "jgap/io/log/CurrentLogger.hpp"
#include "jgap/io/tabgap/TabGapIO.hpp"
#include "jgap/jgap.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"

using namespace jgap;
using namespace jgap::utils;

int main(int argc, char** argv) {
    CurrentLogger::initDefault({});
    const auto start = std::chrono::steady_clock::now();

    const std::string xml_filename = (argc > 1) ? argv[1] : "gap_xml_files/gp_WTaCrV.xml";
    const std::string train_xyz = (argc > 2) ? argv[2] : "db_WTaCrV.xyz";
    const std::string output_prefix = (argc > 3) ? argv[3] : "WTaCrV";
    const double ram_limit_gb = 2.0;

    JGAP_LOG_INFO("Loading QUIP XML template from: {}", xml_filename);
    ValuePtr<Potential> pot = QuipXmlConverter::transform(xml_filename);
    auto* gap_ptr = dynamic_cast<GapPotential*>(pot.get());
    if (!gap_ptr) {
        JGAP_LOG_AND_THROW("Parsed XML potential '{}' is not a GapPotential", xml_filename);
    }
    GapPotential& potential = *gap_ptr;

    JGAP_LOG_INFO("Loading training trajectory from: {}", train_xyz);
    const std::vector<Atoms> training_data = Atoms::readAtoms(train_xyz);
    JGAP_LOG_INFO("Loaded {} training frames", training_data.size());

    // ===== Calculate dimensions and estimate RAM =====
    size_t M = 0;
    for (const auto& comp: potential.getComponents()) {
        M += comp->nSparsePoints();
    }

    size_t N = 0;
    for (const auto& frame: training_data) {
        if (frame.getEnergy().has_value()) N += 1;
        if (frame.getForces().has_value()) N += 3 * frame.nAtoms();
        if (frame.getVirials().has_value()) N += 6;
    }

    const double bytes_per_row = static_cast<double>(M) * sizeof(Real);
    const double kb_per_row = bytes_per_row / 1024.0;
    const double mm_bytes = static_cast<double>(M) * bytes_per_row;
    const double mm_mb = mm_bytes / (1024.0 * 1024.0);
    const double mm_gb = mm_bytes / (1024.0 * 1024.0 * 1024.0);
    const double nm_bytes = static_cast<double>(N + M) * bytes_per_row;
    const double nm_mb = nm_bytes / (1024.0 * 1024.0);
    const double nm_gb = nm_bytes / (1024.0 * 1024.0 * 1024.0);
    const double max_bytes = ram_limit_gb * 1024.0 * 1024.0 * 1024.0;

    // In-place peak streaming memory breakdown:
    // 1. workspace: (M + B) * M * 8 bytes (sizeof(Real))
    // 2. A_chunk input buffer: B * M * 8 bytes
    // Total Peak RAM = (M + B) * M * 8 bytes + B * M * 8 bytes
    //                = mm_bytes + 2 * B * bytes_per_row
    size_t max_chunk_rows = 1;
    if (max_bytes > mm_bytes + 2.0 * bytes_per_row) {
        max_chunk_rows = std::max<size_t>(1, static_cast<size_t>((max_bytes - mm_bytes) / (2.0 * bytes_per_row)));
    }
    const double chunk_mb = static_cast<double>(max_chunk_rows) * bytes_per_row / (1024.0 * 1024.0);
    const double peak_streaming_bytes =
        static_cast<double>(M + max_chunk_rows) * bytes_per_row + static_cast<double>(max_chunk_rows) * bytes_per_row;
    const double peak_streaming_gb = peak_streaming_bytes / (1024.0 * 1024.0 * 1024.0);

    std::cout << "\n================ Memory Estimation ================\n";
    std::cout << "Total GAP components            : " << potential.getComponents().size() << "\n";
    std::cout << "Total columns (sparse points) M : " << M << "\n";
    std::cout << "Total training rows N           : " << N << "\n";
    std::cout << "RAM per row                     : " << bytes_per_row << " B (" << std::fixed << std::setprecision(2)
              << kb_per_row << " KB)\n";
    std::cout << "Full (N+M)xM matrix (QRGapFit)  : " << nm_mb << " MB (" << std::setprecision(3) << nm_gb << " GB)\n";
    std::cout << "MxM covariance accumulator      : " << mm_mb << " MB (" << std::setprecision(4) << mm_gb << " GB)\n";
    std::cout << "Minimum in-place QR RAM         : " << std::setprecision(3) << mm_gb << " GB (1 * MxM)\n";
    std::cout << "User RAM limit                  : " << std::setprecision(2) << ram_limit_gb << " GB\n";
    std::cout << "Max chunk capacity (B rows)     : " << max_chunk_rows << " rows (" << std::setprecision(2) << chunk_mb
              << " MB)\n";
    std::cout << "Estimated Peak Streaming RAM    : " << std::setprecision(3) << peak_streaming_gb << " GB\n";
    std::cout << "===================================================\n\n";

    // ===== setup regularization from QUIP command line parameters =====
    const PerConfigTypeRegularizationRules regularization(
        PerConfigTypeSigmas(0.002, 0.1, 0.2),
        "isolated_atom:0.0001:0.04:0.01:0.0:liquid:0.01:0.5:2.0:0.0:"
        "liquid_composition:0.01:0.5:2.0:0.0:liquid_quaternary:0.01:0.5:2.0:0.0:"
        "surf_liquid:0.01:0.4:0.2:0.0:dimer:0.1:1.0:1.0:0.0:short_range:0.05:0.8:0.8:0.0:"
        "composition:0.005:0.1:0.5:0.0:binary_alloys:0.005:0.1:0.5:0.0:"
        "alloy_ints:0.005:0.1:0.5:0.0:quaternary_composition:0.005:0.1:0.5:0.0:"
        "vac_saddle_alloy:0.005:0.1:0.5:0.0:alloy_vacs:0.005:0.1:0.5:0.0:"
        "alloy_edge:0.005:0.1:0.5:0.0:mcmd:0.005:0.1:0.5:0.0:mcmd_ints:0.005:0.1:0.5:0.0:"
        "mcmd_vacs:0.005:0.1:0.5:0.0:alloy_surf:0.01:0.4:0.5:0.0:"
        "bcc_distorted_dense:0.005:0.1:0.5:0.0:fcc:0.01:0.1:0.5:0.0:hcp:0.01:0.1:0.5:0.0:"
        "C15:0.01:0.1:0.5:0.0:sc:0.01:0.1:0.5:0.0:dia:0.01:0.1:0.5:0.0:"
        "surface_100:0.005:0.1:0.5:0.0:surface_110:0.005:0.1:0.5:0.0:surface_111:0.005:0.1:0.5:0.0:"
        "surface_112:0.005:0.1:0.5:0.0:damian_binaries:0.005:0.1:0.5:0.0:binary_A15:0.01:0.1:1.0:0.0:"
        "binary_C14:0.01:0.1:1.0:0.0:binary_C15:0.01:0.1:1.0:0.0:binary_C36:0.01:0.1:1.0:0.0"
    );

    // ===== fit using Streaming QR (or Split QR) =====
    JGAP_LOG_INFO("Fitting HEA GAP potential using StreamingQrGapFit with limit {} GB...", ram_limit_gb);
    // StreamingQrGapFit fitter(1e-8, ram_limit_gb);

    // Elemental QR with automatic single-species splits
    ElementalQRGapFit fitter(1e-8, ram_limit_gb);

    auto sigmas = regularization.determineForAll(training_data);
    fitter.fit(potential, training_data, sigmas);

    const std::string potential_file = output_prefix + ".jgap.h5";
    SerializationRegistry<Potential>::serialize(pot, potential_file);
    JGAP_LOG_INFO("Saved fitted potential to {}", potential_file);

    // ===== tabulate potential =====
    JGAP_LOG_INFO("Tabulating fitted GAP potential...");
    standardTabulation(potential, output_prefix);

    std::cout << "Execution time: " << formatDuration(elapsedMillisSince(start)) << std::endl;
    return 0;
}
