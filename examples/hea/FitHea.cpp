#include <chrono>
#include <iostream>
#include <string>
#include <vector>

#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/core/atomic/io/XYZData.hpp"
#include "jgap/core/fit/gap/regularization/PerConfigTypeRegularizationRules.hpp"
#include "jgap/core/potentials/gap/GapPotential.hpp"
#include "jgap/core/potentials/tabgap/TabGapPotential.hpp"
#include "jgap/experimental/fit/gap/StreamingQrGapFit.hpp"
#include "jgap/io/convert/QuipXmlConverter.hpp"
#include "jgap/io/log/CurrentLogger.hpp"
#include "jgap/io/tabgap/TabGapIO.hpp"
#include "jgap/jgap.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"

using namespace jgap;

int main(int argc, char** argv) {
    CurrentLogger::initDefault({});
    const auto start = std::chrono::steady_clock::now();

    const std::string xml_filename = (argc > 1) ? argv[1] : "gap_xml_files/gp_WTaCrV.xml";
    const std::string train_xyz = (argc > 2) ? argv[2] : "db_WTaCrV.xyz";
    const std::string output_prefix = (argc > 3) ? argv[3] : "WTaCrV";

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

    // ===== setup regularization from QUIP command line parameters =====
    const PerConfigTypeRegularizationRules regularization(
        ConfigSigmas(0.002, 0.1, 0.2),
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

    // ===== fit using Streaming QR =====
    JGAP_LOG_INFO("Fitting HEA GAP potential using StreamingQrGapFit...");
    StreamingQrGapFit fitter(1e-8, 500);
    fitter.fit(potential, training_data, regularization);

    const std::string potential_file = output_prefix + ".jgap.h5";
    SerializationRegistry<Potential>::serialize(pot, potential_file);
    JGAP_LOG_INFO("Saved fitted potential to {}", potential_file);

    // ===== tabulate potential =====
    JGAP_LOG_INFO("Tabulating fitted GAP potential...");
    TabulationData tabulation_data = potential.tabulate(
        {.max_cutoffs = potential.getCutoffs(), .max_eam_density = 10.0, .n_grid_2b = 5000, .n_grid_3b = {80, 80, 80}}
    );
    TabGapPotential tabgap{tabulation_data};

    const Filenames tabgap_files = TabGapIO::write(tabgap, output_prefix);
    JGAP_LOG_INFO("Saved tabGAP to: {}", vectorToString(tabgap_files));

    std::cout << "Execution time: " << formatDuration(elapsedMillisSince(start)) << std::endl;
    return 0;
}
