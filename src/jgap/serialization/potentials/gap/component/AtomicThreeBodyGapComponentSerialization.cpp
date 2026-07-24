#include "AtomicThreeBodyGapComponentSerialization.hpp"

namespace jgap {

    using AtomicThreeBodySerialization1 = AtomicThreeBodyGapComponentSerialization<2, SquaredExpKernel<1, 1>>;
    using AtomicThreeBodySerialization2 = AtomicThreeBodyGapComponentSerialization<3, SquaredExpKernel<2, 1>>;
    using AtomicThreeBodySerialization3 = AtomicThreeBodyGapComponentSerialization<4, SquaredExpKernel<3, 1>>;
    using AtomicThreeBodySerialization4 = AtomicThreeBodyGapComponentSerialization<5, SquaredExpKernel<4, 1>>;

    template class AtomicThreeBodyGapComponentSerialization<2, SquaredExpKernel<1, 1>>;
    template class AtomicThreeBodyGapComponentSerialization<3, SquaredExpKernel<2, 1>>;
    template class AtomicThreeBodyGapComponentSerialization<4, SquaredExpKernel<3, 1>>;
    template class AtomicThreeBodyGapComponentSerialization<5, SquaredExpKernel<4, 1>>;

    REGISTER_SERIALIZATION(AtomicThreeBodySerialization1, GapComponent);
    REGISTER_SERIALIZATION(AtomicThreeBodySerialization2, GapComponent);
    REGISTER_SERIALIZATION(AtomicThreeBodySerialization3, GapComponent);
    REGISTER_SERIALIZATION(AtomicThreeBodySerialization4, GapComponent);
}
