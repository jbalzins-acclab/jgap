#ifndef JGAP_TABULATABLE_HPP
#define JGAP_TABULATABLE_HPP

#include "data/tabulation/TabulationData.hpp"

namespace jgap {
    class Tabulatable {
    public:
        virtual ~Tabulatable() = default;
        virtual void tabulate(TabulationData &table) = 0;
    };
}

#endif