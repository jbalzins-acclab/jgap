#include "core/potentials/CompositePotential.hpp"

#include "io/log/CurrentLogger.hpp"

namespace jgap {
    std::shared_ptr<CompositePotential> CompositePotential::fromDataNode(const DataNode &params) {

        JGAP_LOG_DEBUG("Parsing composite potential");

        std::map<std::string, std::shared_ptr<Potential>> potentials{};
        for (auto &[label, potentialParams] : REQUIRE(params, "potentials").asObject()) {
            potentials[label] = REGISTRY_GET(Potential, potentialParams);
        }

        return std::make_shared<CompositePotential>(potentials);
    }

    CompositePotential::CompositePotential(const std::map<std::string, std::shared_ptr<Potential>> &potentialsMap)
        : potentials(potentialsMap) {
    }

    DataNode CompositePotential::serialize() {
        DataNode result{};
        result["potentials"] = DataNode::object();
        for (const auto &[label, potential] : potentials) {
            result["potentials"][label] = potential->serialize();
            result["potentials"][label]["type"] = potential->getType();
        }
        return result;
    }

    CutoffRanges CompositePotential::getCutoff() {
        CutoffRanges res{};
        for (const auto &potential: potentials | std::views::values) {
            res += potential->getCutoff();
        }
        return res;
    }

    Predictions CompositePotential::predict(const AtomicStructure &structure) {
        Predictions result{};
        for (const auto &potential : potentials | std::views::values) {
            result = result + potential->predict(structure);
        }
        return result;
    }

    void CompositePotential::tabulate(TabulationData &table) {
        for (const auto &[label, potential] : potentials) {
            JGAP_LOG_DEBUG("Tabulating {} potential", label);
            potential->tabulate(table);
        }
    }
}
