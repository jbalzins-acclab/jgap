#ifndef JGAP_CALCULATIONTYPE_HPP
#define JGAP_CALCULATIONTYPE_HPP

namespace jgap {
    /**
     * Indicator if force/virial(or other derivative) calculations are needed,
     * intended to be used as a template variable to avoid repetition of otherwise similar calculations.
     *
     * For example, when @link NeighbourList finds instances of @link Cluster in a structure
     * @link CalculationType#ValueOnly will prevent it from calculating derivatives irrelevant for e.g., sparsification,
     * meanwhile, @link CalculationType#WithGradients will ensure forces and virials are calculated when needed.
     */
    enum class CalculationType {
        ValueOnly,
        WithGradients
    };

    using CalculationType::ValueOnly;
    using CalculationType::WithGradients;
}

#endif
