#include "ManyBodyGapComponentSerialization.hpp"

namespace jgap {

    // EAM-style many-body component: ManyBodyGapComponent<1, SquaredExpKernel<1, 0>>
    using ManyBody1Serialization = ManyBodyGapComponentSerialization<1, SquaredExpKernel<1, 0>>;

    template class ManyBodyGapComponentSerialization<1, SquaredExpKernel<1, 0>>;

    REGISTER_SERIALIZATION(ManyBody1Serialization, GapComponent);
}
