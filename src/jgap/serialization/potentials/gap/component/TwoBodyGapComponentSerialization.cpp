#include "TwoBodyGapComponentSerialization.hpp"

#include "jgap/experimental/kernels/CauchyKernel.hpp"
#include "jgap/experimental/kernels/WendlandKernel.hpp"

namespace jgap {

    using TwoBodySerialization1 = TwoBodyGapComponentSerialization<2, SquaredExpKernel<1, 1>>;
    using TwoBodySerialization2 = TwoBodyGapComponentSerialization<3, SquaredExpKernel<2, 1>>;
    using TwoBodySerialization3 = TwoBodyGapComponentSerialization<4, SquaredExpKernel<3, 1>>;

    using TwoBodyCauchySerialization1 = TwoBodyGapComponentSerialization<2, CauchyKernel<1, 1>>;
    using TwoBodyCauchySerialization2 = TwoBodyGapComponentSerialization<3, CauchyKernel<2, 1>>;
    using TwoBodyCauchySerialization3 = TwoBodyGapComponentSerialization<4, CauchyKernel<3, 1>>;

    using TwoBodyWendlandSerialization1 = TwoBodyGapComponentSerialization<2, WendlandKernel<1, 1>>;
    using TwoBodyWendlandSerialization2 = TwoBodyGapComponentSerialization<3, WendlandKernel<2, 1>>;
    using TwoBodyWendlandSerialization3 = TwoBodyGapComponentSerialization<4, WendlandKernel<3, 1>>;

    template class TwoBodyGapComponentSerialization<2, SquaredExpKernel<1, 1>>;
    template class TwoBodyGapComponentSerialization<3, SquaredExpKernel<2, 1>>;
    template class TwoBodyGapComponentSerialization<4, SquaredExpKernel<3, 1>>;

    template class TwoBodyGapComponentSerialization<2, CauchyKernel<1, 1>>;
    template class TwoBodyGapComponentSerialization<3, CauchyKernel<2, 1>>;
    template class TwoBodyGapComponentSerialization<4, CauchyKernel<3, 1>>;

    template class TwoBodyGapComponentSerialization<2, WendlandKernel<1, 1>>;
    template class TwoBodyGapComponentSerialization<3, WendlandKernel<2, 1>>;
    template class TwoBodyGapComponentSerialization<4, WendlandKernel<3, 1>>;

    REGISTER_SERIALIZATION(TwoBodySerialization1, GapComponent);
    REGISTER_SERIALIZATION(TwoBodySerialization2, GapComponent);
    REGISTER_SERIALIZATION(TwoBodySerialization3, GapComponent);

    REGISTER_SERIALIZATION(TwoBodyCauchySerialization1, GapComponent);
    REGISTER_SERIALIZATION(TwoBodyCauchySerialization2, GapComponent);
    REGISTER_SERIALIZATION(TwoBodyCauchySerialization3, GapComponent);

    REGISTER_SERIALIZATION(TwoBodyWendlandSerialization1, GapComponent);
    REGISTER_SERIALIZATION(TwoBodyWendlandSerialization2, GapComponent);
    REGISTER_SERIALIZATION(TwoBodyWendlandSerialization3, GapComponent);
}
