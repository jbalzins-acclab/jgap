//
// Created by Jegors Balzins on 9.7.2025.
//
#include "core/potentials/CompositePotential.hpp"

namespace jgap {
    CompositePotential::CompositePotential(map<string, shared_ptr<Potential>> &potentials) {
        _potentials = std::move(potentials);
    }

    CompositePotential::CompositePotential(const nlohmann::json &params) {
        _potentials = {};
        for (const auto &[label, potentialParams] : params["potentials"].items()) {
            _potentials[label] = ParserRegistry<Potential>::get(potentialParams);
        }
    }

    nlohmann::json CompositePotential::serialize() {
        nlohmann::json result{};
        result["potentials"] = nlohmann::json{};
        for (const auto &[label, potential] : _potentials) {
            result["potentials"][label] = potential->serialize();
            result["potentials"][label]["type"] = potential->getType();
        }
        return result;
    }

    double CompositePotential::getCutoff() {
        double cutoff = 0.0;
        for (const auto &potential: _potentials | views::values) {
            cutoff = max(cutoff, potential->getCutoff());
        }
        return cutoff;
    }

    PotentialPrediction CompositePotential::predict(const AtomicStructure &structure) {
        PotentialPrediction result{};
        for (const auto &potential : _potentials | views::values) {
            // auto pp =  potential->predict(structure);
            result = result + potential->predict(structure);
        }
        return result;
    }
}
