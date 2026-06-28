// Example: read a serialized potential and predict on a structure file.
//
// Usage: read_and_predict <pot.h5> <in.xyz> <out.xyz>
//   deserializes <pot.h5> (any potential type), predicts energy/forces/virials for every frame in
//   <in.xyz>, and writes the annotated frames to <out.xyz>.

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "core/atomic/Atoms.hpp"
#include "core/potentials/Potential.hpp"
#include "core/ValuePtr.hpp"
#include "serialization/SerializationRegistry.hpp"
#include "utils/Utils.hpp"
#include "io/log/CurrentLogger.hpp"

using namespace jgap;

int main(int argc, char** argv) {
    CurrentLogger::initDefault({.stdout_log_debug = true});

    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <pot.h5> <in.xyz> <out.xyz>\n";
        return 1;
    }
    const std::string potential_file = argv[1];
    const std::string in_xyz = argv[2];
    const std::string out_xyz = argv[3];

    // Deserialize without knowing the concrete type: the registry picks the right deserializer.
    const ValuePtr<Potential> potential = SerializationRegistry<Potential>::deserialize(potential_file);

    std::vector<Atoms> frames = Atoms::readAtoms(in_xyz);

    std::ofstream out(out_xyz);
    if (!out.is_open()) {
        JGAP_LOG_AND_THROW("Could not open {} for writing", out_xyz);
    }
    for (Atoms& atoms : frames) {
        atoms << potential->calculateEnergy(atoms);
        atoms.write(out);
    }

    JGAP_LOG_INFO("Wrote predictions for {} frame(s) to {}", frames.size(), out_xyz);
    return 0;
}
