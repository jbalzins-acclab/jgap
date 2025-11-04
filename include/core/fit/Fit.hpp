#ifndef JGAP_FIT_HPP
#define JGAP_FIT_HPP

#include "data/Vector3.hpp"
#include "core/potentials/Potential.hpp"

#include <memory>
#include <core/tabulation/Tabulation.hpp>

using namespace std;

namespace jgap {
    class Fit {
    public:
        virtual ~Fit() = default;
        Fit(const nlohmann::json& params) {

            if (!params.contains("tabulation")) return;

            nlohmann::json tabulationParams = params["tabulation"];
            if (!tabulationParams.contains("type")) deduceTabulationType(tabulationParams);

            optionalTabulation = ParserRegistry<Tabulation>::get(tabulationParams);
        }

        shared_ptr<Potential> fit(const vector<AtomicStructure>& trainingData) {
            auto potential = fitWithoutTabulation(trainingData);

            if (optionalTabulation == nullptr) return potential;

            return optionalTabulation->tabulate(potential);
        }

    protected:
        virtual shared_ptr<Potential> fitWithoutTabulation(const vector<AtomicStructure>& trainingData) = 0;

    private:
        shared_ptr<Tabulation> optionalTabulation = nullptr;

        void deduceTabulationType(nlohmann::json& tabulationParams) {

        }
    };
}

#endif
