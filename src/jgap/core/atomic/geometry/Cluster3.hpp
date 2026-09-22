#ifndef JGAP_CLUSTER3_HPP
#define JGAP_CLUSTER3_HPP

#include <concepts>
#include <cstddef>
#include "jgap/core/atomic/geometry/Separation.hpp"

namespace jgap {
    struct Cluster3 {
        size_t idx0{};
        size_t idx1{};
        size_t idx2{};
        Separation separation01{};
        Separation separation02{};
        Separation separation12{};
    };

    template<typename Func>
    concept Cluster3Callback = std::invocable<Func, const Cluster3&>;
}

#endif
