//
// Created by Jegors Balzins on 19.6.2025.
//

#ifndef EAMDENSITYCALCULATOR_HPP
#define EAMDENSITYCALCULATOR_HPP

#include "data/BasicDataTypes.hpp"
#include "data/kernels/EamKernelIndex.hpp"
#include "data/params/EamDescriptorParams.hpp"
#include "pair_functions/EamPairFunction.hpp"

namespace jgap {
    class EamDensityCalculator {
    public:
        explicit EamDensityCalculator(EamDescriptorParams params);
        ~EamDensityCalculator() = default;

        [[nodiscard]] EamKernelIndex calculate(const AtomicStructure &structure) const;

        [[nodiscard]] double getCutoff() const { return _cutoff; }

    private:
        double _cutoff;
        shared_ptr<EamPairFunction> _defaultPairFunction;
        map<OrderedSpeciesPair, shared_ptr<EamPairFunction>> _pairFunctions;

        static shared_ptr<EamPairFunction> selectPairFunction(EamDensityCalculationParams params,
                                                              optional<OrderedSpeciesPair> orderedSpeciesPair);
    };
}

#endif //EAMDENSITYCALCULATOR_HPP
