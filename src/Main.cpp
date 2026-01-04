#include <exception>
#include <iostream>
#include <ostream>

#include "core/descriptors/3b/ThreeBodyDescriptorFinder.hpp"
#include "core/kernels/ThreeBodyKernelCollection.hpp"
#include "core/sparsification/HistogramUniformSparsifier.hpp"
#include "utils/Utils.hpp"

// #include "app/CLIApp.hpp"

using namespace jgap;
using namespace std;

int main(int argc, char** argv) {
#ifndef DEBUG
    try {
        //return jgap::CLIApp::main(argc, argv);

        JGAP_LOG_WARN("TIMESTEP1");
        auto data = readXyz("resources/xyz-samples/feni-train.xyz", 5.0);

        JGAP_LOG_WARN("TIMESTEP2");
        std::shared_ptr<ThreeBodyDescriptorFinder> finder = std::make_shared<ThreeBodyDescriptorFinder>(4.0, 0.6, false);

        map<EncodedSpeciesSets, std::vector<Descriptor<3, 3>>> pts{};
        for (const auto& as: data) {
            auto mp = finder->find(as);
            for ( const auto& [k,v]: mp) {
                if (!pts.contains(k)) pts[k] = {};
                auto &arr = pts[k];
                for (auto& h: v) arr.push_back(h);
            }
        }

        JGAP_LOG_WARN("TIMESTEP3");
        std::map<EncodedSpeciesSets, std::vector<std::shared_ptr<RealKernel<3, 3>>>> kernels_per_species_triplet;
        auto a = HistogramUniformSparsifier<3, 3>(30, 500);
        for (auto& [sp, ppts]: pts) {
            auto sparse_pts = a.selectSparsePoints(ppts);
            kernels_per_species_triplet[sp] = {};
            auto& b = kernels_per_species_triplet[sp];
            for (auto& pt: sparse_pts) {
                std::shared_ptr<RealKernel<3, 3>> kernel = std::make_shared<SquaredExpKernel<3, 3>>(1.0, array{1.0, 1.0, 1.0}, pt.value, pt.f_cut);
                b.push_back(kernel);
            }
        }
        finder->calculate_derivatives = true;

        JGAP_LOG_WARN("TIMESTEP4");
        ThreeBodyKernelCollection collection(finder, kernels_per_species_triplet);
        long long n = 0;
        double r = 0;
        for (const auto& as: data) {
            auto aaa = collection.covariate(as).back();
            n += aaa.forces_optional.size();
            r += aaa.forces_optional.back().len();
        }
        JGAP_LOG_WARN("TIMESTEP5: {}, {}", n, r);

        return EXIT_SUCCESS;

    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }
#else
    // Not to lose the exception when doing a debug run.
    return jgap::CLIApp::main(argc, argv);
#endif
}