#include <exception>
#include <iostream>
#include <ostream>
#include <oneapi/tbb/parallel_for_each.h>

#include "core/descriptors/3b/ThreeBodyDescriptorFinder.hpp"
#include "core/kernels/ThreeBodyKernelCollection.hpp"
#include "core/sparsification/HistogramUniformSparsifier.hpp"
#include "utils/Utils.hpp"

// #include "app/CLIApp.hpp"

using namespace jgap;
using namespace std;

int main(int argc, char** argv) {
//#ifndef DEBUG
    try {
        //return jgap::CLIApp::main(argc, argv);

        JGAP_LOG_WARN("TIMESTEP1");
        auto data = readXyz("resources/xyz-samples/feni-train.xyz", 4.0);

        

        /*JGAP_LOG_WARN("TIMESTEP2");
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
        auto a = make_shared<HistogramUniformSparsifier<3, 3>>(30, 500);
        for (auto& [sp, ppts]: pts) {
            auto sparse_pts = a->selectSparsePoints(ppts);
            kernels_per_species_triplet[sp] = {};
            auto& b = kernels_per_species_triplet[sp];
            JGAP_LOG_INFO("encoded: {} {} {} {} {} {}::{}", sp[0], sp[1], sp[2], sp[3], sp[4], sp[5], ppts.size());
            for (auto& pt: sparse_pts) {
                std::shared_ptr<RealKernel<3, 3>> kernel = std::make_shared<SquaredExpKernel<3, 3>>(1.0, array{1.0, 1.0, 1.0}, pt.value, pt.f_cut);
                b.push_back(kernel);
                //JGAP_LOG_INFO("{}:{}:{}", pt.value[0], pt.value[1], pt.value[2]);
            }
        }
        finder->calculate_derivatives = true;

        JGAP_LOG_WARN("TIMESTEP4");
        ThreeBodyKernelCollection collection(finder, kernels_per_species_triplet);
        atomic_int64_t n = 0;
        double r = 0;
        atomic_int64_t i = 0;

        tbb::parallel_for_each(data.begin(), data.end(), [&](AtomicStructure& as)-> void {
        //for (auto& as: data) {
                //JGAP_LOG_WARN("iii{}", i.load());
            if (i++ % 100 == 0) {
                JGAP_LOG_WARN("iii{}", i.load());
            }
            auto aaa = collection.covariate(as).back();
            n += aaa.forces_optional.size();
            if (!aaa.forces_optional.empty()) {
                r += aaa.forces_optional.back().len();
            }
        }
        );
        JGAP_LOG_WARN("TIMESTEP5: {}, {}", n.load(), r);

        return EXIT_SUCCESS;
        */

    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }
/*#else
    // Not to lose the exception when doing a debug run.
    return jgap::CLIApp::main(argc, argv);
#endif*/
}