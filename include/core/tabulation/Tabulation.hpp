#ifndef JGAP_TABGAPFIT_HPP
#define JGAP_TABGAPFIT_HPP

#include <memory>
#include <nlohmann/json.hpp>

#include "core/potentials/Potential.hpp"
#include "core/potentials/TabGapPotential.hpp"

using namespace std;

namespace jgap {
    class Tabulation {
    public:
        static constexpr string TYPE = "default";

		Tabulation(const nlohmann::json& params);

        shared_ptr<TabGapPotential> tabulate(const shared_ptr<Potential>& potential,
                                             const optional<string>& outputFileNamePrefix);

    private:
        TabulationParams tabulationParams;
    };
}

#endif
