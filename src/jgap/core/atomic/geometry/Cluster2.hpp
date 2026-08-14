#ifndef JGAP_CLUSTER2_HPP
#define JGAP_CLUSTER2_HPP

#include <cstddef>
#include "jgap/core/atomic/geometry/Separation.hpp"

namespace jgap {
    struct Cluster2 {
        size_t idx0{};
        size_t idx1{};
        Separation separation01{};
    };

    template<typename Func>
    concept Cluster2Callback = std::invocable<Func, const Cluster2&>;
}

#endif
