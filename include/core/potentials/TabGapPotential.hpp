//
// Created by balzinsj on 17/10/2025.
//

#ifndef JGAP_TABGAPPOTENTIAL_HPP
#define JGAP_TABGAPPOTENTIAL_HPP

namespace jgap {
    class TabGapPotential {
    public:
        TabGapPotential(nlohmann::json& params);
    };

    REGISTER_PARSER("tabgap", Potential, TabGapPotential)
}

#endif // JGAP_TABGAPPOTENTIAL_HPP
