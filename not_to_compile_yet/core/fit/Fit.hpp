#ifndef JGAP_FIT_HPP
#define JGAP_FIT_HPP

#include "data/Vector3.hpp"
#include "core/potentials/Potential.hpp"

#include <memory>
#include <core/tabulation/SimpleTabulation.hpp>

namespace jgap {

    using DataSet = std::vector<AtomicStructure>;

    class Fit {
    public:
        virtual ~Fit() = default;
        virtual std::shared_ptr<Potential> fit(const std::vector<AtomicStructure>& training_data) = 0;
    };
}

#endif
