#ifndef JGAP_EAMPAIRFUNCTION_HPP
#define JGAP_EAMPAIRFUNCTION_HPP

#include <map>
#include <set>
#include "../TwoBodyTransformation.hpp"
#include "jgap/core/ValuePtr.hpp"
#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/core/atomic/neighbours/NeighbourLists.hpp"
#include "jgap/core/potentials/gap/component/ManyBodyGapComponent.hpp"

namespace jgap {
    enum class EamMode { FSsym, FSgen, EAM, Blind };

    class EamPairFunction : public TwoBodyTransformation<1> {
    public:
        Cutoffs getCutoffs() const override { return Cutoffs{{2, cutoff}}; }

        void setPrefactor(Real p) { prefactor = p; }

        Real getPrefactor() const { return prefactor; }

        Real getCutoff() const { return cutoff; }

        EamPairFunction* clone() const override = 0;

    protected:
        Real cutoff;
        Real prefactor;

        EamPairFunction(Real cutoff = 0.0, Real prefactor = 1.0) : cutoff(cutoff), prefactor(prefactor) {}
    };

}

#endif
