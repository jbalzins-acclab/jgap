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
#include "core/cutoff/CosCutoff.hpp"
#include "core/cutoff/PerriotPolynomialCutoff.hpp"
#include "core/kernels/SquaredExpKernel.hpp"
#include "core/potentials/gap/component/GapComponent.hpp"
#include "core/potentials/gap/GapPotential.hpp"
#include "core/sparsification/HistogramUniformSparsifier.hpp"
#include "core/transform/3b/Angle3bTransformation.hpp"
#include "core/fit/gap/GapFit.hpp"
#include "core/fit/gap/QRGapFit.hpp"
#include "core/fit/gap/regularization/PerConfigTypeRegularizationRules.hpp"
#include "core/fit/gap/regularization/SimpleRegularizationRules.hpp"
#include "core/potentials/CompositePotential.hpp"
#include "core/potentials/gap/component/NBodyGapComponent.hpp"
#include "core/potentials/gap/component/ManyBodyGapComponent.hpp"
#include "core/transform/2b/TwoBodyTransformation.hpp"
#include "core/transform/aggregated/TransformationAggregator.hpp"
#include "core/transform/aggregated/TransformationAggregatorImpl.hpp"
#include "core/transform/eam/FSGenPairFunction.hpp"
#include "utils/Utils.hpp"
#include "utils/convert/QuipXmlConverter.hpp"

using namespace jgap;
using namespace std;

int main(int argc, char** argv) {
    //tbb::global_control control(tbb::global_control::max_allowed_parallelism, 1);

    auto start_time = std::chrono::high_resolution_clock::now();

    /*
    pugi::xml_document quipDocument;
    quipDocument.load_file("gap.xml");
    //quipDocument.load_file("/Users/jegorsbalzins/jgap/notes/etal/sample_baseline.xml");

    auto pot = QuipXmlConverter::transform(quipDocument.document_element()).release();

    auto all = readAtoms("db_train.xyz");
    auto to_be_pred = readAtoms("db_train.xyz")[2];
    auto pred1 = pot->calculateEnergy(to_be_pred);

    //std::println("{} {}");
    //std::cout << pred1.value << " "

    auto* gap = dynamic_cast<GapPotential*>(pot);

    auto reg = PerConfigTypeRegularizationRules(
            0.002,
            0.1,
            0.2,
            0.2,
            {{"isolated_atom", 0.001}},
            {});

    QRGapFit().fit(*gap, readAtoms("db_train.xyz"), reg);
    auto pred2 = gap->calculateEnergy(to_be_pred);

    auto x = dynamic_cast<ManyBodyGapComponent<1>*>(gap->components[0].release());
    auto y = dynamic_cast<SquaredExpKernel<1, 0>*>(x->kernel.release());
    auto z = dynamic_cast<TransformationAggregatorImpl<1, 2>*>(x->aggregator.release());

    return 0;*/
    JGAP_LOG_INFO("Start");
    auto training_data = readAtoms("/Users/jegorsbalzins/jgap/resources/xyz-samples/db_Fe.xyz");
    //auto training_data = readAtoms("/Users/jegorsbalzins/jgap/resources/xyz-samples/feni-train.xyz");

    // ====================================================================================
    // ManyBodyGapComponent with FSGenPairFunction
    // ====================================================================================
    auto eam_pf = FSGenPairFunction(4.5, 3.0);
    auto neighbour_lists_eam = NeighbourList::generateFor(training_data, eam_pf);

    std::vector<Descriptor<1>> descs_eam{};
    for (auto &nl: neighbour_lists_eam) {
        std::map<size_t, Descriptor<1>> agg_descs;
        auto raw = nl.findAllClusters(SpeciesSet<2, HasCentralAtom>("Fe", "Fe"));
        for (auto &cluster: raw) {
            auto val = eam_pf.evaluate(cluster);
            agg_descs[cluster.atom_indexes[0]].value[0] += val.value[0];
        }
        for (const auto& [idx, desc] : agg_descs) {
            descs_eam.push_back(desc);
        }
    }

    auto sparsifier_eam = HistogramUniformSparsifier<1>(120, 20);
    auto sp_eam = sparsifier_eam.selectSparsePoints(descs_eam);

    auto kernel_eam = SquaredExpKernel<1, 0>(1.0, {1.0});

    auto eam_aggregator = std::make_unique<TransformationAggregatorImpl<1, 2>>("Fe");
    eam_aggregator->extend(SpeciesSet<2, HasCentralAtom>("Fe", "Fe"), std::make_unique<FSGenPairFunction>(eam_pf));

    auto component_eam = std::make_unique<ManyBodyGapComponent<1>>(
        std::move(eam_aggregator),
        std::make_unique<SquaredExpKernel<1, 0>>(kernel_eam),
        sp_eam
    );

    // ====================================================================================
    // 3-Body and 2-Body Components
    // ====================================================================================
    auto cutoff3 = CosCutoff(3.7, 0.6);
    auto cutoff2 = CosCutoff(4.5, 1.0);

    auto trans3 = Angle3bTransformation(std::make_unique<CosCutoff>(cutoff3));
    auto trans2 = TwoBodyTransformation(std::make_unique<CosCutoff>(cutoff2));

    auto kernel3 = SquaredExpKernel<3, 1>(1.0, {1.0, 1.0, 1.0});
    auto kernel2 = SquaredExpKernel<1, 1>(10.0, {1.0});

    auto sparsifier3 = HistogramUniformSparsifier<3+1>(100, 500);
    auto sparsifier2 = HistogramUniformSparsifier<1+1>(120, 20);

    auto neighbour_lists3 = NeighbourList::generateFor(training_data, trans3);
    std::vector<Descriptor<4>> descs{};
    for (auto &nl: neighbour_lists3) {
        auto raw = nl.findAllClusters(SpeciesSet<3, HasCentralAtom>("Fe", "Fe", "Fe"));
        for (auto &cluster: raw) {
            descs.push_back(trans3.evaluate(cluster));
        }
    }
    auto sp3 = sparsifier3.selectSparsePoints(descs);

    auto neighbour_lists2 = NeighbourList::generateFor(training_data, trans2);
    std::vector<Descriptor<2>> descs2{};
    for (auto &nl: neighbour_lists2) {
        auto raw = nl.findAllClusters(SpeciesSet<2, Symmetric>("Fe", "Fe"));
        for (auto &cluster: raw) {
            descs2.push_back(trans2.evaluate(cluster));
        }
    }
    auto sp2 = sparsifier2.selectSparsePoints(descs2);

    auto component3 = std::make_unique<NBodyGapComponent<4, 3, HasCentralAtom, SquaredExpKernel<3, 1>>>(
            SpeciesSet<3, HasCentralAtom>("Fe", "Fe", "Fe"),
            std::make_unique<Angle3bTransformation>(std::move(trans3)),
            //std::make_unique<SquaredExpKernel<3, 1>>(kernel3),
            kernel3,
            sp3
        );

    auto component2 = std::make_unique<NBodyGapComponent<2, 2, Symmetric, SquaredExpKernel<1, 1>>>(
            SpeciesSet<2, Symmetric>("Fe", "Fe"),
            std::make_unique<TwoBodyTransformation>(std::move(trans2)),
            //std::make_unique<SquaredExpKernel<1, 1>>(kernel2),
            kernel2,
            sp2
        );

    // ====================================================================================
    // Construct Potential
    // ====================================================================================
    std::vector<GapComponent::Ptr> components;
    // Add all components
    components.push_back(std::move(component_eam));
    components.push_back(std::move(component2));
    components.push_back(std::move(component3));

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