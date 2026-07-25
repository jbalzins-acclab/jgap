#include "ManyBodyGapComponentSerialization.hpp"

#include "jgap/experimental/kernels/CauchyKernel.hpp"
#include "jgap/experimental/kernels/WendlandKernel.hpp"

namespace jgap {

    // EAM-style many-body component: ManyBodyGapComponent<1, SquaredExpKernel<1, 0>>
    using ManyBody1Serialization = ManyBodyGapComponentSerialization<1, SquaredExpKernel<1, 0>>;
    using ManyBody2Serialization = ManyBodyGapComponentSerialization<2, SquaredExpKernel<2, 0>>;
    using ManyBody3Serialization = ManyBodyGapComponentSerialization<3, SquaredExpKernel<3, 0>>;
    using ManyBody4Serialization = ManyBodyGapComponentSerialization<4, SquaredExpKernel<4, 0>>;

    using ManyBodyCauchy1Serialization = ManyBodyGapComponentSerialization<1, CauchyKernel<1, 0>>;
    using ManyBodyCauchy2Serialization = ManyBodyGapComponentSerialization<2, CauchyKernel<2, 0>>;
    using ManyBodyCauchy3Serialization = ManyBodyGapComponentSerialization<3, CauchyKernel<3, 0>>;
    using ManyBodyCauchy4Serialization = ManyBodyGapComponentSerialization<4, CauchyKernel<4, 0>>;

    using ManyBodyWendland1Serialization = ManyBodyGapComponentSerialization<1, WendlandKernel<1, 0>>;
    using ManyBodyWendland2Serialization = ManyBodyGapComponentSerialization<2, WendlandKernel<2, 0>>;
    using ManyBodyWendland3Serialization = ManyBodyGapComponentSerialization<3, WendlandKernel<3, 0>>;
    using ManyBodyWendland4Serialization = ManyBodyGapComponentSerialization<4, WendlandKernel<4, 0>>;

    template class ManyBodyGapComponentSerialization<1, SquaredExpKernel<1, 0>>;
    template class ManyBodyGapComponentSerialization<2, SquaredExpKernel<2, 0>>;
    template class ManyBodyGapComponentSerialization<3, SquaredExpKernel<3, 0>>;
    template class ManyBodyGapComponentSerialization<4, SquaredExpKernel<4, 0>>;

    template class ManyBodyGapComponentSerialization<1, CauchyKernel<1, 0>>;
    template class ManyBodyGapComponentSerialization<2, CauchyKernel<2, 0>>;
    template class ManyBodyGapComponentSerialization<3, CauchyKernel<3, 0>>;
    template class ManyBodyGapComponentSerialization<4, CauchyKernel<4, 0>>;

    template class ManyBodyGapComponentSerialization<1, WendlandKernel<1, 0>>;
    template class ManyBodyGapComponentSerialization<2, WendlandKernel<2, 0>>;
    template class ManyBodyGapComponentSerialization<3, WendlandKernel<3, 0>>;
    template class ManyBodyGapComponentSerialization<4, WendlandKernel<4, 0>>;

    REGISTER_SERIALIZATION(ManyBody1Serialization, GapComponent);
    REGISTER_SERIALIZATION(ManyBody2Serialization, GapComponent);
    REGISTER_SERIALIZATION(ManyBody3Serialization, GapComponent);
    REGISTER_SERIALIZATION(ManyBody4Serialization, GapComponent);

    REGISTER_SERIALIZATION(ManyBodyCauchy1Serialization, GapComponent);
    REGISTER_SERIALIZATION(ManyBodyCauchy2Serialization, GapComponent);
    REGISTER_SERIALIZATION(ManyBodyCauchy3Serialization, GapComponent);
    REGISTER_SERIALIZATION(ManyBodyCauchy4Serialization, GapComponent);

    REGISTER_SERIALIZATION(ManyBodyWendland1Serialization, GapComponent);
    REGISTER_SERIALIZATION(ManyBodyWendland2Serialization, GapComponent);
    REGISTER_SERIALIZATION(ManyBodyWendland3Serialization, GapComponent);
    REGISTER_SERIALIZATION(ManyBodyWendland4Serialization, GapComponent);
}
