#ifndef JGAP_REGULARIZATION_HPP
#define JGAP_REGULARIZATION_HPP
#include <optional>
#include <vector>

#include "jgap/core/Real.hpp"
#include "jgap/core/atomic/energy/Virials.hpp"

#include "../../../Vector3.hpp"

namespace jgap {
    struct Regularization {
        std::optional<Real> energy;
        std::optional<Virials> virials;
        std::optional<std::vector<Vector3>> forces;
    };
}

#endif
