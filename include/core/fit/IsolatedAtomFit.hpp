#ifndef JGAP_ISOLATEDATOMFIT_HPP
#define JGAP_ISOLATEDATOMFIT_HPP

#include "Fit.hpp"

namespace jgap {
    class IsolatedAtomFit : public Fit {
    public:
        static constexpr string TYPE = "isolated_atom";

        IsolatedAtomFit(const nlohmann::json& params);
        shared_ptr<Potential> fit(const vector<AtomicStructure> &trainingData) override;

    private:
        bool _errorOnUnknownSpecies;
    };

    REGISTER_PARSER(Fit, IsolatedAtomFit)
}

#endif
