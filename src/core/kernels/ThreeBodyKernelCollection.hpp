#ifndef JGAP_THREEBODYKERNELCOLLECTION_HPP
#define JGAP_THREEBODYKERNELCOLLECTION_HPP
#include <ranges>

#include "KernelCollection.hpp"
#include "core/cutoff/CutoffFunction.hpp"
#include "core/descriptors/3b/ThreeBodyDescriptorFinder.hpp"

namespace jgap {

    class ThreeBodyKernelCollection : public KernelCollection/*, Tabulatable, Serializable*/ {
    public:


        std::vector<Predictions> covariate(const Box &atomic_structure) override {

        }

        std::vector<std::shared_ptr<Matrix>> selfCovariate() override {
        }

        std::vector<std::shared_ptr<IKernel>> getKernels() override {
        }


    private:
    };
}

#endif
