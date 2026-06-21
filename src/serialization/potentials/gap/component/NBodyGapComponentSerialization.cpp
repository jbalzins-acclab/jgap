#include "NBodyGapComponentSerialization.hpp"

namespace jgap {

    // x-numbers+cutoff <- 2-body
    using NBodySerialization12 = NBodyGapComponentSerialization<2, 2, FullSymmetry, SquaredExpKernel<1, 1>>;
    using NBodySerialization22 = NBodyGapComponentSerialization<3, 2, FullSymmetry, SquaredExpKernel<2, 1>>;
    using NBodySerialization32 = NBodyGapComponentSerialization<4, 2, FullSymmetry, SquaredExpKernel<3, 1>>;

    // x-numbers+cutoff <- 3-body
    using NBodySerialization13 = NBodyGapComponentSerialization<2, 3, HasCentralAtom, SquaredExpKernel<1, 1>>;
    using NBodySerialization23 = NBodyGapComponentSerialization<3, 3, HasCentralAtom, SquaredExpKernel<2, 1>>;
    using NBodySerialization33 = NBodyGapComponentSerialization<4, 3, HasCentralAtom, SquaredExpKernel<3, 1>>;
    using NBodySerialization43 = NBodyGapComponentSerialization<5, 3, HasCentralAtom, SquaredExpKernel<4, 1>>;

    template class NBodyGapComponentSerialization<2, 2, FullSymmetry, SquaredExpKernel<1, 1>>;
    template class NBodyGapComponentSerialization<3, 2, FullSymmetry, SquaredExpKernel<2, 1>>;
    template class NBodyGapComponentSerialization<4, 2, FullSymmetry, SquaredExpKernel<3, 1>>;

    template class NBodyGapComponentSerialization<2, 3, HasCentralAtom, SquaredExpKernel<1, 1>>;
    template class NBodyGapComponentSerialization<3, 3, HasCentralAtom, SquaredExpKernel<2, 1>>;
    template class NBodyGapComponentSerialization<4, 3, HasCentralAtom, SquaredExpKernel<3, 1>>;
    template class NBodyGapComponentSerialization<5, 3, HasCentralAtom, SquaredExpKernel<4, 1>>;

    REGISTER_SERIALIZATION(NBodySerialization12, GapComponent);
    REGISTER_SERIALIZATION(NBodySerialization22, GapComponent);
    REGISTER_SERIALIZATION(NBodySerialization32, GapComponent);
    REGISTER_SERIALIZATION(NBodySerialization13, GapComponent);
    REGISTER_SERIALIZATION(NBodySerialization23, GapComponent);
    REGISTER_SERIALIZATION(NBodySerialization33, GapComponent);
    REGISTER_SERIALIZATION(NBodySerialization43, GapComponent);
}
