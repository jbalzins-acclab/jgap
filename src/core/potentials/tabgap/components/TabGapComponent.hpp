#ifndef JGAP_TABGAPCOMPONENT_HPP
#define JGAP_TABGAPCOMPONENT_HPP
#include "core/atomic/energy/AtomicQuantity.hpp"
#include "core/atomic/neighbours/NeighbourList.hpp"
#include "core/tabulation/TabulationData.hpp"

namespace jgap {
    class TabGapComponent {
    public:
        virtual ~TabGapComponent() = default;
        virtual AtomicQuantity energy(const NeighbourList& nl) const = 0;
        virtual Cutoffs getCutoffs() const = 0;
        virtual void tabulate(TabulationData& tables) const = 0;
        virtual std::set<Species> getAllSpecies() const = 0;

        virtual TabGapComponent* clone() const = 0;
    };
}

#endif
