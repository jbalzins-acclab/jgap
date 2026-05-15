#ifndef JGAP_POTENTIAL_HPP
#define JGAP_POTENTIAL_HPP

#include <string>
#include <functional>

#include "core/atomic/energy/AtomicQuantity.hpp"
#include "data/Cutoffs.hpp"
#include "io/Serializable.hpp"

namespace jgap {
    class Potential {
    public:
        virtual ~Potential() = default;

        virtual AtomicQuantity energy(const Atoms &atoms) = 0;

        virtual Cutoffs getCutoffs() = 0;
        //virtual void tabulate(TabulationData& table) = 0;
    };
}

#endif
