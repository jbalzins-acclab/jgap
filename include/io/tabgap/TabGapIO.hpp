#ifndef JGAP_TABGAPIO_HPP
#define JGAP_TABGAPIO_HPP

#include "core/potentials/TabGapPotential.hpp"

namespace jgap {
    class TabGapIO {
    public:
        static void read(shared_ptr<TabGapPotential> potential);
        static void write(shared_ptr<TabGapPotential> potential);

    private:
        static void writeH5(shared_ptr<TabGapPotential> potential);
        static void writeEamFs(shared_ptr<TabGapPotential> potential, size_t index);

        static void readH5(shared_ptr<TabGapPotential> potential);
        static void readEamFs(shared_ptr<TabGapPotential> potential);
    };
}

#endif