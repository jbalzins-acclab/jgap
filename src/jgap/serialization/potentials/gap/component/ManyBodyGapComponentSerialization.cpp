#include "ManyBodyGapComponentSerialization.hpp"

namespace jgap {

    // EAM-style many-body component: ManyBodyGapComponent<1, SquaredExpKernel<1, 0>>
    using ManyBody1Serialization = ManyBodyGapComponentSerialization<1, SquaredExpKernel<1, 0>>;
    using ManyBody2Serialization = ManyBodyGapComponentSerialization<2, SquaredExpKernel<2, 0>>;
    using ManyBody3Serialization = ManyBodyGapComponentSerialization<3, SquaredExpKernel<3, 0>>;
    using ManyBody4Serialization = ManyBodyGapComponentSerialization<4, SquaredExpKernel<4, 0>>;

    template class ManyBodyGapComponentSerialization<1, SquaredExpKernel<1, 0>>;
    template class ManyBodyGapComponentSerialization<2, SquaredExpKernel<2, 0>>;
    template class ManyBodyGapComponentSerialization<3, SquaredExpKernel<3, 0>>;
    template class ManyBodyGapComponentSerialization<4, SquaredExpKernel<4, 0>>;

    REGISTER_SERIALIZATION(ManyBody1Serialization, GapComponent);
    REGISTER_SERIALIZATION(ManyBody2Serialization, GapComponent);
    REGISTER_SERIALIZATION(ManyBody3Serialization, GapComponent);
    REGISTER_SERIALIZATION(ManyBody4Serialization, GapComponent);
}
