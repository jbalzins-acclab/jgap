#include "core/potentials/CompositePotential.hpp"

#include "io/log/CurrentLogger.hpp"

namespace jgap {
    CompositePotential::CompositePotential(const map<string, shared_ptr<Potential>>& potentials) {
        _potentials = potentials;
    }

    CompositePotential::CompositePotential(const nlohmann::json &params) {
        JGAP_LOG_DEBUG("Parsing composite potential");
        _potentials = {};
        for (const auto &[label, potentialParams] : require(params, "potentials").items()) {
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

    CutoffRanges CompositePotential::getCutoff() {
        CutoffRanges res{};
        for (const auto &potential: _potentials | views::values) {
            res += potential->getCutoff();
        }
        return res;
    }

    Predictions CompositePotential::predict(const AtomicStructure &structure) {
        Predictions result{};
        for (const auto &potential : _potentials | views::values) {
            result = result + potential->predict(structure);
        }
        return result;
    }

    void CompositePotential::tabulate(TabulationData &table) {
        for (const auto &[label, potential] : _potentials) {
            JGAP_LOG_DEBUG("Tabulating {} potential", label);
            potential->tabulate(table);
        }
    }
}
