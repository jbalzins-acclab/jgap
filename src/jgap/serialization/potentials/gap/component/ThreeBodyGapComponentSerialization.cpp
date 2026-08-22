#include "ThreeBodyGapComponentSerialization.hpp"

#include "jgap/experimental/kernels/CauchyKernel.hpp"
#include "jgap/experimental/kernels/WendlandKernel.hpp"

namespace jgap {

    using ThreeBodySerialization1 = ThreeBodyGapComponentSerialization<2, SquaredExpKernel<1, 1>>;
    using ThreeBodySerialization2 = ThreeBodyGapComponentSerialization<3, SquaredExpKernel<2, 1>>;
    using ThreeBodySerialization3 = ThreeBodyGapComponentSerialization<4, SquaredExpKernel<3, 1>>;
    using ThreeBodySerialization4 = ThreeBodyGapComponentSerialization<5, SquaredExpKernel<4, 1>>;

    using ThreeBodyCauchySerialization1 = ThreeBodyGapComponentSerialization<2, CauchyKernel<1, 1>>;
    using ThreeBodyCauchySerialization2 = ThreeBodyGapComponentSerialization<3, CauchyKernel<2, 1>>;
    using ThreeBodyCauchySerialization3 = ThreeBodyGapComponentSerialization<4, CauchyKernel<3, 1>>;
    using ThreeBodyCauchySerialization4 = ThreeBodyGapComponentSerialization<5, CauchyKernel<4, 1>>;

    using ThreeBodyWendlandSerialization1 = ThreeBodyGapComponentSerialization<2, WendlandKernel<1, 1>>;
    using ThreeBodyWendlandSerialization2 = ThreeBodyGapComponentSerialization<3, WendlandKernel<2, 1>>;
    using ThreeBodyWendlandSerialization3 = ThreeBodyGapComponentSerialization<4, WendlandKernel<3, 1>>;
    using ThreeBodyWendlandSerialization4 = ThreeBodyGapComponentSerialization<5, WendlandKernel<4, 1>>;

    template class ThreeBodyGapComponentSerialization<2, SquaredExpKernel<1, 1>>;
    template class ThreeBodyGapComponentSerialization<3, SquaredExpKernel<2, 1>>;
    template class ThreeBodyGapComponentSerialization<4, SquaredExpKernel<3, 1>>;
    template class ThreeBodyGapComponentSerialization<5, SquaredExpKernel<4, 1>>;

    template class ThreeBodyGapComponentSerialization<2, CauchyKernel<1, 1>>;
    template class ThreeBodyGapComponentSerialization<3, CauchyKernel<2, 1>>;
    template class ThreeBodyGapComponentSerialization<4, CauchyKernel<3, 1>>;
    template class ThreeBodyGapComponentSerialization<5, CauchyKernel<4, 1>>;

    template class ThreeBodyGapComponentSerialization<2, WendlandKernel<1, 1>>;
    template class ThreeBodyGapComponentSerialization<3, WendlandKernel<2, 1>>;
    template class ThreeBodyGapComponentSerialization<4, WendlandKernel<3, 1>>;
    template class ThreeBodyGapComponentSerialization<5, WendlandKernel<4, 1>>;

    REGISTER_SERIALIZATION(ThreeBodySerialization1, GapComponent);
    REGISTER_SERIALIZATION(ThreeBodySerialization2, GapComponent);
    REGISTER_SERIALIZATION(ThreeBodySerialization3, GapComponent);
    REGISTER_SERIALIZATION(ThreeBodySerialization4, GapComponent);

    REGISTER_SERIALIZATION(ThreeBodyCauchySerialization1, GapComponent);
    REGISTER_SERIALIZATION(ThreeBodyCauchySerialization2, GapComponent);
    REGISTER_SERIALIZATION(ThreeBodyCauchySerialization3, GapComponent);
    REGISTER_SERIALIZATION(ThreeBodyCauchySerialization4, GapComponent);

    REGISTER_SERIALIZATION(ThreeBodyWendlandSerialization1, GapComponent);
    REGISTER_SERIALIZATION(ThreeBodyWendlandSerialization2, GapComponent);
    REGISTER_SERIALIZATION(ThreeBodyWendlandSerialization3, GapComponent);
    REGISTER_SERIALIZATION(ThreeBodyWendlandSerialization4, GapComponent);
}
