#include "core/potentials/CompositePotential.hpp"

#include "io/log/CurrentLogger.hpp"

namespace jgap {
    CompositePotential::CompositePotential(const map<string, shared_ptr<Potential>>& potentials) {
        _potentials = potentials;
    }

    CompositePotential::CompositePotential(const nlohmann::json &params) {
        CurrentLogger::get()->debug("Parsing composite potential");
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

    Predictions CompositePotential::predict(const AtomicStructure &structure) {
        Predictions result{};
        for (const auto &potential : _potentials | views::values) {
            // auto pp =  potential->predict(structure);
            result = result + potential->predict(structure);
        }
        return result;
    }

    void CompositePotential::tabulate(TabulationData &table) {
        for (const auto &[label, potential] : _potentials) {
            CurrentLogger::get()->debug("Tabulating {} potential", label);
            potential->tabulate(table);
        }
    }
}
