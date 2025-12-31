#ifndef JGAP_POTENTIALWITHMETA_HPP
#define JGAP_POTENTIALWITHMETA_HPP

#include <memory>

#include "core/potentials/Potential.hpp"

namespace jgap {
    struct PotentialWithMeta {
        std::shared_ptr<Potential> potential;
        bool saved = false;
        bool to_be_saved = false;
    };
}

#endif