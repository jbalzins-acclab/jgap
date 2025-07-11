#ifndef ISOLATEDATOMFIT_HPP
#define ISOLATEDATOMFIT_HPP

#include "Fit.hpp"

namespace jgap {
    class IsolatedAtomFit : public Fit {
    public:
        ~IsolatedAtomFit() override = default;

        explicit IsolatedAtomFit(const nlohmann::json& params);
        string getType() override { return "isolated_atom"; }

        shared_ptr<Potential> fit(const vector<AtomicStructure> &trainingData) override;

    private:
        bool _errorOnUnknownSpecies;
    };

    REGISTER_PARSER("isolated_atom", Fit, IsolatedAtomFit)
}

#endif //ISOLATEDATOMFIT_HPP
