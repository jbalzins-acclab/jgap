#ifndef THREEBODYKERNELINDEX_HPP
#define THREEBODYKERNELINDEX_HPP

#include <vector>

using namespace std;

namespace jgap {
    struct ThreeBodyKernelIndexEntity {
        array<size_t, 3> atomIndex;

        array<Vector3, 3> r_ij;
        array<Vector3, 3> grad_rij_wrt_rj;
        // TODO: add 3rd if needed
        double fCut01;
        double fCut02;
        double dfCut_dr_01;
        double dfCut_dr_02;

        Vector3 q;
        array<Vector3, 3> dq_k_dr_ij;
    };

    using ThreeBodyKernelIndex = vector<ThreeBodyKernelIndexEntity>;

    struct ThreeBodyDescriptorData {
        Vector3 q;
        double fCut;
    };
}

#endif //THREEBODYKERNELINDEX_HPP
