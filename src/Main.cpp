#include <exception>
#include <iostream>
#include <ostream>
#include <ranges>
#include <memory>
#include <oneapi/tbb/parallel_for_each.h>
#include <oneapi/tbb/global_control.h>

#include "core/atomic/io/XYZData.hpp"
#include "core/atomic/Atoms.hpp"
#include "core/cutoff/CosCutoff.hpp"
#include "core/cutoff/PerriotPolynomialCutoff.hpp"
#include "core/kernels/SquaredExpKernel.hpp"
#include "core/kernels/groups/BlindKernelGroup.hpp"
#include "core/kernels/groups/PerSpeciesKernelGroup.hpp"
#include "core/potentials/gap/GapComponent.hpp"
#include "core/potentials/gap/GapPotential.hpp"
#include "core/sparsification/HistogramUniformSparsifier.hpp"
#include "core/transform/3b/Angle3bTransformer.hpp"
#include "core/fit/gap/GapFit.hpp"
#include "core/fit/gap/QRGapFit.hpp"
#include "core/fit/gap/regularization/SimpleRegularizationRules.hpp"
#include "utils/Utils.hpp"


// #include "app/CLIApp.hpp"

using namespace jgap;
using namespace std;

int main(int argc, char** argv) {
    /*tbb::global_control control(
        tbb::global_control::max_allowed_parallelism, 1
    );*/

    auto training_data = readAtoms("/Users/jegorsbalzins/jgap/resources/xyz-samples/db_Fe.xyz");

    auto cutoff1 = CosCutoff(3.7, 0.6);
    auto cutoff2 = PerriotPolynomialCutoff(3.5, 3.7);

    std::shared_ptr<Transformer<4, 3>> trans1 = std::make_shared<Angle3bTransformer<CosCutoff>>(cutoff1);
    std::shared_ptr<Transformer<4, 3>> trans2 = std::make_shared<Angle3bTransformer<CosCutoff>>(cutoff1);

    auto kernel1 = SquaredExpKernel<3, 1, 3>(3.0, {1.0, 1.0, 1.0});
    auto kernel2 = SquaredExpKernel<3, 1, 3>(2.0, {0.5, 0.5, 0.5});

    auto sparsifier1 = std::make_shared<HistogramUniformSparsifier<4>>(100, 500);
    auto sparsifier2 = std::make_shared<HistogramUniformSparsifier<4>>(120, 500);

    // Generate neighbour lists using the helper
    auto neighbour_lists = NeighbourList::generateFor(training_data, trans1, trans2);

    // Construct Kernel Groups
    std::shared_ptr<KernelGroup<4, 3>> kernel_group1 = std::make_shared<BlindKernelGroup<SquaredExpKernel<3, 1, 3>>>(
        BlindKernelGroup<SquaredExpKernel<3, 1, 3>>::constructFromTrainingData(
            kernel1, trans1, sparsifier1, neighbour_lists
        )
    );

    std::map<SpeciesSet, std::shared_ptr<Sparsifier<4>>> empty_sparsifiers;
    std::map<SpeciesSet, SquaredExpKernel<3, 1, 3>> empty_kernels;

    auto kernel_group2 = std::make_shared<PerSpeciesKernelGroup<SquaredExpKernel<3, 1, 3>>>(
        PerSpeciesKernelGroup<SquaredExpKernel<3, 1, 3>>::constructFromTrainingData(
            trans2,
            neighbour_lists,
            empty_sparsifiers,
            sparsifier2,
            empty_kernels,
            kernel2,
            true
        )
    );

    // Construct Gap Components

    auto component1 = std::make_shared<GapComponent<4, 3>>(trans1, kernel_group1);
    auto component2 = std::make_shared<GapComponent<4, 3>>(trans2, kernel_group2);

    std::vector<std::shared_ptr<IGapComponent>> components = {component2};

    // Construct Potential
    auto potential = std::make_unique<GapPotential>(std::move(components));

    std::shared_ptr<RegularizationRules> regularization_rules = std::make_shared<SimpleRegularizationRules>();

    // Fit
    QRGapFit fitter;
    fitter.fit(*potential, training_data, neighbour_lists, {}, regularization_rules);

}