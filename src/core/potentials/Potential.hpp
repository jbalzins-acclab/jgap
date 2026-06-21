#ifndef JGAP_POTENTIAL_HPP
#define JGAP_POTENTIAL_HPP

#include <string>
#include <functional>

#include "core/atomic/energy/AtomicQuantity.hpp"
#include "core/atomic/Atoms.hpp"
#include "Cutoffs.hpp"
#include "core/tabulation/TabulationData.hpp"

namespace jgap {
    class Potential {
    public:
        virtual ~Potential() = default;

        virtual AtomicQuantity calculateEnergy(const Atoms &atoms) const = 0;

        virtual Cutoffs getCutoffs() const = 0;

        virtual Potential* clone() const = 0;

        virtual void fillTables(TabulationData& table) const = 0;

        virtual TabulationData tabulate(TabulationParams tabulation_params) const {
            TabulationData tabulation_data{tabulation_params};

            fillTables(tabulation_data);

            return tabulation_data;
        }
    };
}

#endif
