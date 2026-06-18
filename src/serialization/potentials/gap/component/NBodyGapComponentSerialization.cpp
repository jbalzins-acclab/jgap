#include "NBodyGapComponentSerialization.hpp"

namespace jgap {

    // 2-body: NBodyGapComponent<2, 2, Symmetric, SquaredExpKernel<1, 1>>
    using NBody2BodySerialization = NBodyGapComponentSerialization<2, 2, Symmetric, SquaredExpKernel<1, 1>>;
    // 3-body: NBodyGapComponent<4, 3, HasCentralAtom, SquaredExpKernel<3, 1>>
    using NBody3BodySerialization = NBodyGapComponentSerialization<4, 3, HasCentralAtom, SquaredExpKernel<3, 1>>;

    template class NBodyGapComponentSerialization<2, 2, Symmetric, SquaredExpKernel<1, 1>>;
    template class NBodyGapComponentSerialization<4, 3, HasCentralAtom, SquaredExpKernel<3, 1>>;

    REGISTER_SERIALIZATION(NBody2BodySerialization, GapComponent);
    REGISTER_SERIALIZATION(NBody3BodySerialization, GapComponent);
}
