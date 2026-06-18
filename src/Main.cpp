#include <exception>
#include <iostream>
#include <ostream>
#include <ranges>
#include <memory>
#include <chrono>
#include <iomanip>
#include <pugixml.hpp>
#include <oneapi/tbb/parallel_for_each.h>
#include <oneapi/tbb/global_control.h>

#include "core/atomic/io/XYZData.hpp"
#include "core/atomic/Atoms.hpp"
#include "io/tabgap/TabGapIO.hpp"
#include "utils/Utils.hpp"
#include "utils/convert/QuipXmlConverter.hpp"
#include "utils/gap/StandardGapFit.hpp"
#include "core/potentials/tabgap/TabGapPotential.hpp"
#include "core/transform/eam/FSGenPairFunction.hpp"
#include "core/fit/gap/regularization/SimpleRegularizationRules.hpp"

using namespace jgap;
using namespace std;

int main(int argc, char** argv) {
    //tbb::global_control control(tbb::global_control::max_allowed_parallelism, 1);

    CurrentLogger::initDefault({.stdout_log_debug = true});

    auto start_time = std::chrono::high_resolution_clock::now();

    JGAP_LOG_INFO("Start");
    auto training_data = readAtoms("resources/xyz-samples/db_Fe.xyz");

    StandardGapParams params{120};
    params.eam_pf = FSGenPairFunction(4.5, 3.0);
    params.eam_mode = EamMode::Blind;
    params.regularization_rules = SimpleRegularizationRules();

    auto potential = standardGapFit(training_data, params);

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

    return 0;
}