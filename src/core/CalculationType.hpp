#ifndef JGAP_CALCULATIONTYPE_HPP
#define JGAP_CALCULATIONTYPE_HPP

namespace jgap {
    enum class CalculationType {
        ValueOnly,
        WithGradients
    };

    using CalculationType::ValueOnly;
    using CalculationType::WithGradients;
}

#endif
