#ifndef COMPOSITEFIT_HPP
#define COMPOSITEFIT_HPP

#include "Fit.hpp"
#include "core/neighbours/NeighbourFinder.hpp"

namespace jgap {
    class CompositePotentialFit : public Fit {
    public:
        static constexpr string TYPE = "composite";

        CompositePotentialFit(const nlohmann::json& params);

        shared_ptr<Potential> fit(const vector<AtomicStructure> &trainingData) override;

    private:
        optional<shared_ptr<Potential>> _externalPotential;
        map<string, shared_ptr<Fit>> _fits;
        vector<string> _fitOrder;

        static vector<AtomicStructure> subtractExternalContribution(const vector<AtomicStructure> &originalData,
                                                                    const shared_ptr<Potential> &potential);

    protected:
        shared_ptr<Potential> fitWithoutTabulation(const vector<AtomicStructure> &trainingData) override;
    };

    REGISTER_PARSER(Fit, CompositePotentialFit)
}

#endif //COMPOSITEFIT_HPP
