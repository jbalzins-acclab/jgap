#ifndef JGAP_TRANSFORMATIONTRAITS_HPP
#define JGAP_TRANSFORMATIONTRAITS_HPP

#include <array>
#include <cstddef>
#include "3b/Angle3bTransformer.hpp" // Forward declaration might be better if headers get complex

namespace jgap {

    // Generic trait: by default, assumes all dimensions are used for all separations
    template<typename Transformation>
    struct TransformationTraits {
        static constexpr bool is_specialized = false;
    };

    // Specialization for Angle3bTransformation
    template<typename TCutoff>
    struct TransformationTraits<Angle3bTransformation<TCutoff>> {
        static constexpr bool is_specialized = true;

        // For each separation, list the descriptor dimensions that have non-zero derivatives
        static constexpr std::array<std::array<size_t, 3>, 2> non_zero_dims_r01_r02 = {{
            {0, 1, 3}, // Non-zero dims for r_01
            {0, 1, 3}  // Non-zero dims for r_02
        }};

        static constexpr std::array<size_t, 1> non_zero_dims_r12 = {{2}}; // Non-zero dims for r_12
    };

}

#endif //JGAP_TRANSFORMATIONTRAITS_HPP