#include "AtomicThreeBodyGapComponentSerialization.hpp"

#include "jgap/experimental/kernels/CauchyKernel.hpp"
#include "jgap/experimental/kernels/WendlandKernel.hpp"

namespace jgap {

    using AtomicThreeBodySerialization1 = AtomicThreeBodyGapComponentSerialization<2, SquaredExpKernel<1, 1>>;
    using AtomicThreeBodySerialization2 = AtomicThreeBodyGapComponentSerialization<3, SquaredExpKernel<2, 1>>;
    using AtomicThreeBodySerialization3 = AtomicThreeBodyGapComponentSerialization<4, SquaredExpKernel<3, 1>>;
    using AtomicThreeBodySerialization4 = AtomicThreeBodyGapComponentSerialization<5, SquaredExpKernel<4, 1>>;

    using AtomicThreeBodyCauchySerialization1 = AtomicThreeBodyGapComponentSerialization<2, CauchyKernel<1, 1>>;
    using AtomicThreeBodyCauchySerialization2 = AtomicThreeBodyGapComponentSerialization<3, CauchyKernel<2, 1>>;
    using AtomicThreeBodyCauchySerialization3 = AtomicThreeBodyGapComponentSerialization<4, CauchyKernel<3, 1>>;
    using AtomicThreeBodyCauchySerialization4 = AtomicThreeBodyGapComponentSerialization<5, CauchyKernel<4, 1>>;

    using AtomicThreeBodyWendlandSerialization1 = AtomicThreeBodyGapComponentSerialization<2, WendlandKernel<1, 1>>;
    using AtomicThreeBodyWendlandSerialization2 = AtomicThreeBodyGapComponentSerialization<3, WendlandKernel<2, 1>>;
    using AtomicThreeBodyWendlandSerialization3 = AtomicThreeBodyGapComponentSerialization<4, WendlandKernel<3, 1>>;
    using AtomicThreeBodyWendlandSerialization4 = AtomicThreeBodyGapComponentSerialization<5, WendlandKernel<4, 1>>;

    template class AtomicThreeBodyGapComponentSerialization<2, SquaredExpKernel<1, 1>>;
    template class AtomicThreeBodyGapComponentSerialization<3, SquaredExpKernel<2, 1>>;
    template class AtomicThreeBodyGapComponentSerialization<4, SquaredExpKernel<3, 1>>;
    template class AtomicThreeBodyGapComponentSerialization<5, SquaredExpKernel<4, 1>>;

    template class AtomicThreeBodyGapComponentSerialization<2, CauchyKernel<1, 1>>;
    template class AtomicThreeBodyGapComponentSerialization<3, CauchyKernel<2, 1>>;
    template class AtomicThreeBodyGapComponentSerialization<4, CauchyKernel<3, 1>>;
    template class AtomicThreeBodyGapComponentSerialization<5, CauchyKernel<4, 1>>;

    template class AtomicThreeBodyGapComponentSerialization<2, WendlandKernel<1, 1>>;
    template class AtomicThreeBodyGapComponentSerialization<3, WendlandKernel<2, 1>>;
    template class AtomicThreeBodyGapComponentSerialization<4, WendlandKernel<3, 1>>;
    template class AtomicThreeBodyGapComponentSerialization<5, WendlandKernel<4, 1>>;

    REGISTER_SERIALIZATION(AtomicThreeBodySerialization1, GapComponent);
    REGISTER_SERIALIZATION(AtomicThreeBodySerialization2, GapComponent);
    REGISTER_SERIALIZATION(AtomicThreeBodySerialization3, GapComponent);
    REGISTER_SERIALIZATION(AtomicThreeBodySerialization4, GapComponent);

    REGISTER_SERIALIZATION(AtomicThreeBodyCauchySerialization1, GapComponent);
    REGISTER_SERIALIZATION(AtomicThreeBodyCauchySerialization2, GapComponent);
    REGISTER_SERIALIZATION(AtomicThreeBodyCauchySerialization3, GapComponent);
    REGISTER_SERIALIZATION(AtomicThreeBodyCauchySerialization4, GapComponent);

    REGISTER_SERIALIZATION(AtomicThreeBodyWendlandSerialization1, GapComponent);
    REGISTER_SERIALIZATION(AtomicThreeBodyWendlandSerialization2, GapComponent);
    REGISTER_SERIALIZATION(AtomicThreeBodyWendlandSerialization3, GapComponent);
    REGISTER_SERIALIZATION(AtomicThreeBodyWendlandSerialization4, GapComponent);
}
