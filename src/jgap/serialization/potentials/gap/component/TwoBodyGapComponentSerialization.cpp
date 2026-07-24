#include "TwoBodyGapComponentSerialization.hpp"

namespace jgap {

    using TwoBodySerialization1 = TwoBodyGapComponentSerialization<2, SquaredExpKernel<1, 1>>;
    using TwoBodySerialization2 = TwoBodyGapComponentSerialization<3, SquaredExpKernel<2, 1>>;
    using TwoBodySerialization3 = TwoBodyGapComponentSerialization<4, SquaredExpKernel<3, 1>>;

    template class TwoBodyGapComponentSerialization<2, SquaredExpKernel<1, 1>>;
    template class TwoBodyGapComponentSerialization<3, SquaredExpKernel<2, 1>>;
    template class TwoBodyGapComponentSerialization<4, SquaredExpKernel<3, 1>>;

    REGISTER_SERIALIZATION(TwoBodySerialization1, GapComponent);
    REGISTER_SERIALIZATION(TwoBodySerialization2, GapComponent);
    REGISTER_SERIALIZATION(TwoBodySerialization3, GapComponent);
}
