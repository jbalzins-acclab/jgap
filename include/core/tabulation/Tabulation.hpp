#ifndef JGAP_TABULATION_HPP
#define JGAP_TABULATION_HPP

#include <memory>

#include "Tabulation.hpp"
#include "core/potentials/Potential.hpp"
#include "core/potentials/TabGapPotential.hpp"

namespace jgap {
    class Tabulation {
    public:
        std::shared_ptr<TabGapPotential> tabulate(const std::shared_ptr<Potential>& potential);
    };
}

#endif