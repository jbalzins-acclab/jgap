#include <iostream>
#include <chrono>
#include <iomanip>
#include <string>

#include "core/atomic/io/XYZData.hpp"
#include "core/atomic/Atoms.hpp"
#include "io/tabgap/TabGapIO.hpp"
#include "utils/gap/StandardGapFit.hpp"
#include "core/potentials/tabgap/TabGapPotential.hpp"
#include "core/transform/eam/FSGenPairFunction.hpp"
#include "core/fit/gap/regularization/SimpleRegularizationRules.hpp"
#include "core/potentials/Potential.hpp"
#include "core/ValuePtr.hpp"
#include "serialization/SerializationRegistry.hpp"
#include "io/log/CurrentLogger.hpp"

using namespace jgap;
using namespace std;

namespace {
    const std::string TrainingDataFile = "resources/xyz-samples/db_Fe.xyz";
    const std::string PotentialFile = "fe-potential.jgap.h5";
}

void fit() {
    auto start_time = std::chrono::high_resolution_clock::now();

    JGAP_LOG_INFO("Start");
    auto training_data = readAtoms(TrainingDataFile);

    StandardGapParams params{120};
    params.eam_pf = FSGenPairFunction(4.5, 3.0);
    params.eam_mode = EamMode::Blind;
    params.regularization_rules = SimpleRegularizationRules();

    auto potential = standardGapFit(training_data, params);

    // Serialize the fitted potential. The registry stamps the node with the concrete type, so it can be
    // read back (see readAndTest) without the caller knowing what kind of potential it is.
    SerializationRegistry<Potential>::serialize(ValuePtr<Potential>(potential), PotentialFile);
    JGAP_LOG_INFO("Saved fitted potential to {}", PotentialFile);

    TabulationData tabulation_data = potential.tabulate({
            .max_cutoffs = potential.getCutoffs(),
            .max_eam_density = 10.0,
            .n_grid_2b = 5000,
            .n_grid_3b = {80, 80, 80}
        });

    TabGapPotential tabgap{tabulation_data};

    TabGapIO::write(tabgap, "fe-test-1");

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration);
    duration -= minutes;
    auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration);
    duration -= seconds;
    auto milliseconds = duration;

    std::cout << "Execution time: "
              << std::setfill('0') << std::setw(2) << minutes.count() << ":"
              << std::setfill('0') << std::setw(2) << seconds.count() << ":"
              << std::setfill('0') << std::setw(3) << milliseconds.count() << std::endl;

    Atoms to_be_pred = training_data[800];
    to_be_pred.write("orig.xyz");

    to_be_pred << potential.calculateEnergy(to_be_pred);

    to_be_pred.write("gap.xyz");
    to_be_pred << tabgap.calculateEnergy(to_be_pred);

    to_be_pred.write("tabgap.xyz");
}

void readAndTest() {
    JGAP_LOG_INFO("Reading potential from {} (type deduced from file)", PotentialFile);

    // Deserialize without knowing the concrete type: the registry picks the right deserializer.
    ValuePtr<Potential> potential = SerializationRegistry<Potential>::deserialize(PotentialFile);

    auto structures = readAtoms(TrainingDataFile);
    Atoms structure = structures[800];

    AtomicQuantity prediction = potential->calculateEnergy(structure);

    std::cout << "Loaded potential predicted energy on db_Fe.xyz structure #800 ("
              << structure.nAtoms() << " atoms): " << prediction.value << " eV" << std::endl;

    structure << prediction;
    structure.write("loaded-prediction.xyz");
}

int main(int argc, char** argv) {
    CurrentLogger::initDefault({.stdout_log_debug = true});

    fit();
    readAndTest();

    return 0;
}
