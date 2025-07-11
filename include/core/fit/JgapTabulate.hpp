#ifndef JGAPTABULATE_HPP
#define JGAPTABULATE_HPP

#include <memory>

#include "../potentials/JgapPotential.hpp"
#include "../potentials/JtabGapPotential.hpp"

using namespace std;

namespace jgap {

    class JgapTabulate {
    public:
        JgapTabulate(/*TODO*/);
        ~JgapTabulate();

        shared_ptr<void> tabulate(shared_ptr<JgapPotential> jgapPotential);
    };
}

#endif //JGAPTABULATE_HPP
