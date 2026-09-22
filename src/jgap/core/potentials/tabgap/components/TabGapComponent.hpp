#ifndef JGAP_TABGAPCOMPONENT_HPP
#define JGAP_TABGAPCOMPONENT_HPP
#include "jgap/core/atomic/energy/AtomicQuantity.hpp"
#include "jgap/core/atomic/neighbours/NeighbourLists.hpp"
#include "jgap/core/tabulation/TabulationData.hpp"

namespace jgap {
    class TabGapComponent {
    public:
        virtual ~TabGapComponent() = default;
        virtual AtomicQuantity energy(const NeighbourLists& nl) const = 0;
        virtual Cutoffs getCutoffs() const = 0;
        virtual void tabulate(TabulationData& tables) const = 0;
        virtual std::set<Species> getAllSpecies() const = 0;

        virtual TabGapComponent* clone() const = 0;
    };
}

#endif
