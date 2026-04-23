#ifndef JGAP_ANGLE3BTRANSFORMER_HPP
#define JGAP_ANGLE3BTRANSFORMER_HPP

#include "core/cutoff/CutoffFunction.hpp"
#include "core/descriptors/SeparationTransformer.hpp"


namespace jgap {
    class Angle3bTransformer : public SeparationTransformer::FixedDependencies<3, 3+1> {
    public:
        using FixedDependencies::ToBeCalculatedPerDimension;

        std::shared_ptr<CutoffFunction> cutoff;

        std::array<ToBeCalculatedPerDimension, 3+1> evaluate(const TransformFrom &separations) override;
    };
}

#endif
