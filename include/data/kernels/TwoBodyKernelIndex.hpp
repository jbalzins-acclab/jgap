//
// Created by Jegors Balzins on 18.6.2025.
//

#ifndef TWOBODYKERNELINDEX_HPP
#define TWOBODYKERNELINDEX_HPP

namespace jgap {
    struct TwoBodyKernelIndexEntity {
        size_t atomIndex;
        size_t neighbourListIndex;
    };

    using TwoBodyKernelIndex = vector<TwoBodyKernelIndexEntity>;
}

#endif //TWOBODYKERNELINDEX_HPP
