// jgap command-line entry point.
//
// Usage:
//   jgap --predict  <pot.h5> <in.xyz> <out.xyz>
//       Loads the serialized potential from <pot.h5>, predicts energy/forces/virials for every frame in
//       <in.xyz>, and writes the annotated frames to <out.xyz>.
//
//   jgap --tabulate <pot.h5> [key=value ...]
//       Loads the potential, tabulates it (cutoffs taken from the potential; other TabulationParams
//       default unless overridden by key=value args), and writes <pot>.tabgap.h5 plus the EAM
//       <pot>.eam.fs file(s). Recognised keys: r_min_3b, max_eam_density, n_grid_2b, n_grid_3b (the last
//       as three comma-separated counts, e.g. n_grid_3b=80,80,80).

#include <array>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "jgap/core/UnseqFor.hpp"
#include "jgap/core/ValuePtr.hpp"
#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/core/potentials/Potential.hpp"
#include "jgap/core/potentials/tabgap/TabGapPotential.hpp"
#include "jgap/core/tabulation/TabulationParams.hpp"
#include "jgap/io/PotentialLoader.hpp"
#include "../jgap/core/io/log/CurrentLogger.hpp"
#include "jgap/io/tabgap/TabGapIO.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"
#include "jgap/utils/Utils.hpp"

using namespace jgap;

void printUsage(const char* exe) {
    std::cerr << "Usage:\n"
              << "  " << exe << " --predict  <pot.h5> <in.xyz> <out.xyz>\n"
              << "  " << exe
              << " --tabulate <pot.h5> [r_min_3b=..] [max_eam_density=..]"
                 " [n_grid_2b=..] [n_grid_3b=a,b,c]\n";
}

/// Strips a trailing ".h5" so e.g. "pot.h5" -> "pot" (the prefix TabGapIO::write extends).
std::string potentialPrefix(std::string path) {
    if (path.ends_with(".h5")) {
        path.resize(path.size() - 3);
    }
    return path;
}

/// Splits "key=value" into its parts; value is empty when there is no '='.
std::pair<std::string, std::string> splitKeyValue(const std::string& arg) {
    const auto eq = arg.find('=');
    if (eq == std::string::npos) {
        return {arg, ""};
    }
    return {arg.substr(0, eq), arg.substr(eq + 1)};
}

std::array<size_t, 3> parseTriple(const std::string& value) {
    std::array<size_t, 3> out{};
    size_t start = 0;
    for (size_t i = 0; i < 3; i++) {
        const auto comma = value.find(',', start);
        out[i] = std::stoul(value.substr(start, comma - start));
        if (comma == std::string::npos) {
            if (i != 2) JGAP_LOG_AND_THROW("n_grid_3b needs three comma-separated values: {}", value);
            break;
        }
        start = comma + 1;
    }
    return out;
}



int predict(const std::vector<std::string>& args) {
    if (args.size() != 3) {
        std::cerr << "--predict expects: <pot_file> <in.xyz> <out.xyz>\n";
        return 1;
    }
    const std::string& pot_file = args[0];
    const std::string& in_xyz = args[1];
    const std::string& out_xyz = args[2];

    const ValuePtr<Potential> potential = loadPotential(pot_file);

    std::vector<Atoms> frames = Atoms::readAtoms(in_xyz);

    unseqForEach(frames.begin(), frames.end(), [&](Atoms& atoms) {
        atoms.setEnergyAndDerivatives(potential->calculateEnergy(atoms));
    });

    std::ofstream out(out_xyz);
    if (!out.is_open()) {
        JGAP_LOG_AND_THROW("Could not open {} for writing", out_xyz);
    }
    for (const Atoms& atoms: frames) {
        atoms.write(out);
    }

    JGAP_LOG_INFO("Wrote predictions for {} frame(s) to {}", frames.size(), out_xyz);
    return 0;
}

int tabulate(const std::vector<std::string>& args) {
    if (args.empty()) {
        std::cerr << "--tabulate expects: <pot.h5> [key=value ...]\n";
        return 1;
    }
    const std::string& pot_file = args[0];

    const ValuePtr<Potential> potential = SerializationRegistry<Potential>::deserialize(pot_file);

    // Default tabulation params, with cutoffs taken from the potential; override the rest from args.
    TabulationParams params;
    params.max_cutoffs = potential->getCutoffs();

    for (size_t i = 1; i < args.size(); i++) {
        const auto [key, value] = splitKeyValue(args[i]);
        if (key == "r_min_3b") {
            params.r_min_3b = std::stod(value);
        } else if (key == "max_eam_density") {
            params.max_eam_density = std::stod(value);
        } else if (key == "n_grid_2b") {
            params.n_grid_2b = std::stoul(value);
        } else if (key == "n_grid_3b") {
            params.n_grid_3b = parseTriple(value);
        } else {
            JGAP_LOG_AND_THROW("Unknown tabulation parameter '{}'", key);
        }
    }

    const TabGapPotential tabgap{potential->tabulate(params)};
    TabGapIO::write(tabgap, potentialPrefix(pot_file));
    return 0;
}

int main(int argc, char** argv) {
    CurrentLogger::initDefault({});

    if (argc < 2) {
        printUsage(argv[0]);
        return 1;
    }

    const std::string mode = argv[1];
    const std::vector<std::string> args(argv + 2, argv + argc);

    if (mode == "--predict") {
        return predict(args);
    }
    if (mode == "--tabulate") {
        return tabulate(args);
    }

    printUsage(argv[0]);
    return 1;
}
