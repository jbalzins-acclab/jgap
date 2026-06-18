#include "TransformationAggregatorImplSerialization.hpp"

namespace jgap {
    template class TransformationAggregatorImplSerialization<1, 2>;
    template class TransformationAggregatorImplSerialization<1, 3>;

    using AggregatorSerialization1_2 = TransformationAggregatorImplSerialization<1, 2>;
    using AggregatorSerialization1_3 = TransformationAggregatorImplSerialization<1, 3>;
    using AggregatorBase1 = TransformationAggregator<1>;

    REGISTER_SERIALIZATION(AggregatorSerialization1_2, AggregatorBase1);
    REGISTER_SERIALIZATION(AggregatorSerialization1_3, AggregatorBase1);
}