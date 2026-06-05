#ifndef JGAP_EAMPAIRFUNCTION_HPP
#define JGAP_EAMPAIRFUNCTION_HPP

#include "core/transform/ClusterTransformation.hpp"

namespace jgap {
    class EamPairFunction : public ClusterTransformation<1, 2> {
    public:
        Cutoffs getCutoffs() const override {
            return Cutoffs{ {2, cutoff} };
        }

    protected:
        EamPairFunction(Real cutoff = 0.0, Real prefactor = 1.0)
            : cutoff(cutoff), prefactor(prefactor) {}
        Real cutoff;
        Real prefactor;
    };
}

#endif