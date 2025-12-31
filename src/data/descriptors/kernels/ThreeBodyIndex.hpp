#ifndef THREEBODYKERNELINDEX_HPP
#define THREEBODYKERNELINDEX_HPP

#include <vector>

#include "data/Vector3.hpp"

namespace jgap {
    struct ThreeBodyIndexEntity {
        std::array<size_t, 3> atom_index;

        std::array<Vector3, 3> r_ij;
        std::array<Vector3, 3> grad_rij_wrt_rj;
        // TODO: add 3rd if needed
        double fCut01;
        double fCut02;
        double dfCut_dr_01;
        double dfCut_dr_02;

        Vector3 q;
        std::array<Vector3, 3> dq_k_dr_ij;
    };

    using ThreeBodyIndex = std::vector<ThreeBodyIndexEntity>;

    struct ThreeBodyDescriptorData {
        Vector3 q;
        double f_cut;
    };
}

#endif
