#ifndef JGAP_KERNELSET_HPP
#define JGAP_KERNELSET_HPP

#include <vector>

#include "core/kernels/Kernel.hpp"

namespace jgap {
    struct KernelManager {
        std::vector<std::shared_ptr<IKernel>> kernels;
        std::string label;
    };
}

#endif