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

    const std::string xml_filename = (argc > 1) ? argv[1] : "potential_xml_files/gap_HEA.xml";
    const std::string train_xyz = (argc > 2) ? argv[2] : "db/train.xyz";
    const std::string output_prefix = (argc > 3) ? argv[3] : "hea";

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
        "isolated_atom:0.0001:0.04:0.04:0.0:liquid:0.01:0.5:2.0:0.0:dimer:0.01:0.5:2.0:0.0:"
        "short_range:0.01:0.5:2.0:0.0:liquid_surface_100:0.01:0.5:2.0:0.0:"
        "liquid_surface_110:0.01:0.5:2.0:0.0:liquid_surface_111:0.01:0.5:2.0:0.0:"
        "gamma_surface:0.002:0.08:0.5:0.0:liquid_high:0.02:0.8:5.0:0.0:"
        "HEA_liquid:0.01:0.5:2.0:0.0:HEA_short_range:0.01:0.5:2.0:0.0"
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
