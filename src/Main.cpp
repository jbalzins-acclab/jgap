#include <exception>
#include <iostream>
#include <ostream>
#include <ranges>
#include <memory>
#include <chrono>
#include <iomanip>
#include <oneapi/tbb/parallel_for_each.h>
#include <oneapi/tbb/global_control.h>

#include "core/atomic/io/XYZData.hpp"
#include "core/atomic/Atoms.hpp"
#include "core/cutoff/CosCutoff.hpp"
#include "core/cutoff/PerriotPolynomialCutoff.hpp"
#include "core/kernels/SquaredExpKernel.hpp"
#include "core/potentials/gap/component/GapComponent.hpp"
#include "core/potentials/gap/GapPotential.hpp"
#include "core/sparsification/HistogramUniformSparsifier.hpp"
#include "core/transform/3b/DynamicAngle3bTransformation.hpp"
#include "core/fit/gap/GapFit.hpp"
#include "core/fit/gap/QRGapFit.hpp"
#include "core/fit/gap/regularization/SimpleRegularizationRules.hpp"
#include "core/potentials/gap/component/DynamicNBodyGapComponent.hpp"
#include "core/transform/2b/TwoBodyTransformer.hpp"
#include "core/transform/eam/PolycutoffPairFunction.hpp"
#include "utils/Utils.hpp"

using namespace jgap;
using namespace std;

int main(int argc, char** argv) {
    // tbb::global_control control(tbb::global_control::max_allowed_parallelism, 6);

    auto start_time = std::chrono::high_resolution_clock::now();

    JGAP_LOG_INFO("Start");
    auto training_data = readAtoms("/Users/jegorsbalzins/jgap/resources/xyz-samples/db_Fe.xyz");

    /*
    uint64_t space{};
    auto nls = NeighbourList::form(training_data, 4.2);
    for (auto& nl: nls) {
        auto clusters = nl.findAllClusters<3>({"Fe", "Fe", "Fe"});
        std::cout << clusters.size() / nl.nAtoms() << endl;
        if (clusters.size() > 0) {
            space += clusters.size() * sizeof(clusters[0]);
        }
    }
    std::cout << space / 1024 / 1024 << std::endl;
    return 0;
    */

    auto cutoff3 = CosCutoff(3.7, 0.6);
    //auto cutoff2 = PerriotPolynomialCutoff(3.5, 4.5);
    auto cutoff2 = CosCutoff(4.5, 1.0);
    auto pair_function = PolycutoffPairFunction(5.0, 3.0);

    auto trans3 = DynamicAngle3bTransformation(std::make_unique<CosCutoff>(cutoff3));
    auto trans2 = TwoBodyTransformation(cutoff2);
    //auto trans_eam = EamTransformer(pair_function);

    auto kernel3 = SquaredExpKernel<3, 1>(1.0, {1.0, 1.0, 1.0});
    auto kernel2 = SquaredExpKernel<1, 1>(10.0, {1.0});
    auto kernel_eam = SquaredExpKernel<1, 0>(1.0, {1.0});

    auto sparsifier3 = HistogramUniformSparsifier<3+1>(100, 500);
    auto sparsifier2 = HistogramUniformSparsifier<1+1>(120, 20);

    // Generate neighbour lists using the helper
    auto neighbour_lists3 = NeighbourList::generateFor(training_data, trans3);

    std::vector<Descriptor<4>> descs{};
    for (auto &nl: neighbour_lists3) {
        auto raw = nl.findAllClusters<3>({"Fe", "Fe", "Fe"});

        //nl.iterateOverClusters<3>({"Fe", "Fe", "Fe"}, [&](const Cluster<3> &cluster) {
        for (auto &cluster: raw) {
            descs.push_back(trans3.evaluate(cluster));
        }//);
    }

    auto sp3 = sparsifier3.selectSparsePoints(descs);

    auto neighbour_lists2 = NeighbourList::generateFor(training_data, trans2);

    std::vector<Descriptor<2>> descs2{};
    for (auto &nl: neighbour_lists2) {
        auto raw = nl.findAllClusters<2>({"Fe", "Fe"});
        for (auto &cluster: raw) {
            descs2.push_back(trans2.evaluate(cluster));
        }
    }

    auto sp2 = sparsifier2.selectSparsePoints(descs2);

    // Construct Gap Components

    auto component3 = GapComponent::Ptr(new DynamicNBodyGapComponent<4, 3>(
            {"Fe","Fe", "Fe"},
            std::make_unique<DynamicAngle3bTransformation>(std::move(trans3)),
            std::make_unique<SquaredExpKernel<3, 1>>(kernel3),
            sp3
        ));

    auto component2 = GapComponent::Ptr(new DynamicNBodyGapComponent<2, 2>(
            {"Fe","Fe"},
            std::make_unique<TwoBodyTransformation<CosCutoff>>(trans2),
            std::make_unique<SquaredExpKernel<1, 1>>(kernel2),
            sp2
        ));
    // auto component2 = std::make_shared<GapComponent<2, 2>>(trans2, kernel_group2);

    std::vector<GapComponent::Ptr> components;
    components.push_back(std::move(component2));
    components.push_back(std::move(component3));

    // Construct Potential
    auto potential = GapPotential(std::move(components));

    auto regularization_rules = SimpleRegularizationRules();

    // Fit
    QRGapFit fitter;
    fitter.fit(potential, training_data, regularization_rules);

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

    return 0;
}