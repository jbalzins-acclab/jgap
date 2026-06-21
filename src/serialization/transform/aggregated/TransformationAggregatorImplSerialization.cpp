#include "TransformationAggregatorImplSerialization.hpp"

namespace jgap {
    template class TransformationAggregatorImplSerialization<1, 2>;
    template class TransformationAggregatorImplSerialization<1, 3>;
    template class TransformationAggregatorImplSerialization<2, 2>;
    template class TransformationAggregatorImplSerialization<2, 3>;
    template class TransformationAggregatorImplSerialization<3, 2>;
    template class TransformationAggregatorImplSerialization<3, 3>;
    template class TransformationAggregatorImplSerialization<4, 2>;
    template class TransformationAggregatorImplSerialization<4, 3>;

    using AggregatorSerialization1_2 = TransformationAggregatorImplSerialization<1, 2>;
    using AggregatorSerialization1_3 = TransformationAggregatorImplSerialization<1, 3>;
    using AggregatorBase1 = TransformationAggregator<1>;

    using AggregatorSerialization2_2 = TransformationAggregatorImplSerialization<2, 2>;
    using AggregatorSerialization2_3 = TransformationAggregatorImplSerialization<2, 3>;
    using AggregatorBase2 = TransformationAggregator<2>;

    using AggregatorSerialization3_2 = TransformationAggregatorImplSerialization<3, 2>;
    using AggregatorSerialization3_3 = TransformationAggregatorImplSerialization<3, 3>;
    using AggregatorBase3 = TransformationAggregator<3>;

    using AggregatorSerialization4_2 = TransformationAggregatorImplSerialization<4, 2>;
    using AggregatorSerialization4_3 = TransformationAggregatorImplSerialization<4, 3>;
    using AggregatorBase4 = TransformationAggregator<4>;

    REGISTER_SERIALIZATION(AggregatorSerialization1_2, AggregatorBase1);
    REGISTER_SERIALIZATION(AggregatorSerialization1_3, AggregatorBase1);

    REGISTER_SERIALIZATION(AggregatorSerialization2_2, AggregatorBase2);
    REGISTER_SERIALIZATION(AggregatorSerialization2_3, AggregatorBase2);

    REGISTER_SERIALIZATION(AggregatorSerialization3_2, AggregatorBase3);
    REGISTER_SERIALIZATION(AggregatorSerialization3_3, AggregatorBase3);

    REGISTER_SERIALIZATION(AggregatorSerialization4_2, AggregatorBase4);
    REGISTER_SERIALIZATION(AggregatorSerialization4_3, AggregatorBase4);
}