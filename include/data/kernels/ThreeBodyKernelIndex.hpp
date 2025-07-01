#ifndef THREEBODYKERNELINDEX_HPP
#define THREEBODYKERNELINDEX_HPP

namespace jgap {
    struct ThreeBodyKernelIndexEntity {
        size_t atomIndex;
        size_t neighbourListIndex1;
        size_t neighbourListIndex2;
    };

    using ThreeBodyKernelIndex = vector<ThreeBodyKernelIndexEntity>;
}

#endif //THREEBODYKERNELINDEX_HPP
