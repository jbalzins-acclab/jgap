#ifndef JGAP_ISOLATEDATOMFIT_HPP
#define JGAP_ISOLATEDATOMFIT_HPP

#include "Fit.hpp"

namespace jgap {
    class IsolatedAtomFit : public Fit {
    public:
        SETUP_PARSER(Fit, IsolatedAtomFit, isolated_atom)

        IsolatedAtomFit(bool errorOnUnknownSpecies = false) : error_on_unknown_species_(errorOnUnknownSpecies) {}

        std::shared_ptr<Potential> fit(const std::vector<AtomicStructure> &trainingData) override;

    private:
        bool error_on_unknown_species_;
    };
}

#endif
