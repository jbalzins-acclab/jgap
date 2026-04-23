#include "core/fit/CompositePotentialFit.hpp"

#include <fstream>
#include <oneapi/tbb/parallel_for_each.h>

#include "core/potentials/CompositePotential.hpp"
#include "data/DataNode.hpp"

namespace jgap {
    CompositePotentialFit::CompositePotentialFit(const DataNode &params) {
        if (params.contains("external")) {
            JGAP_LOG_INFO("External potential setup");
            _externalPotential = REGISTRY_GET(Potential, params["external"]);
        } else if (params.contains("external_from_file")) {
            JGAP_LOG_INFO(
                format("External potential setup from {}", params["external_from_file"].dump())
            );

            DataNode externalPotentialParams;
            ifstream externalPotentialFile(params["external_from_file"].get<std::string>());
            if (!externalPotentialFile.is_open()) {
                JGAP_LOG_ERROR("Could not open external potential file", true);
            }

            externalPotentialFile >> externalPotentialParams;
            _externalPotential = REGISTRY_GET(Potential, externalPotentialParams);
        } else {
            _externalPotential = {};
        }

        if () {

        }

        _fits = {};
        for (auto &[label, fitParams]: params["fits"].items()) {
            JGAP_LOG_INFO("Picking fitting logic for " + label);
            fitParams[]
            _fits[label] = REGISTRY_GET(Fit, fitParams);
        }

        if (params.contains("fit_order")) {
            for (const auto &label: params["fit_order"]) {
                _fitOrder.push_back(label);
            }
            if (_fitOrder.size() != _fits.size()) {
                JGAP_LOG_WARN("Fit order size mismatches number of fit params");
            }
        } else {
            if (_fits.size() == 1) {
                _fitOrder = {_fits.begin()->first};
            } else if (_fits.size() == 2 && _fits.contains("isolated_atom")) {
                _fitOrder = std::vector(std::views::keys(_fits).begin(), std::views::keys(_fits).end());
                if (_fitOrder[0] != "isolated_atom") {
                    swap(_fitOrder[0], _fitOrder[1]);
                }
            } else {
                JGAP_LOG_ERROR("Fitting order not specified", true);
            }
        }
    }

    std::shared_ptr<Potential> CompositePotentialFit::fit(const std::vector<AtomicStructure> &trainingData) {
        std::vector<AtomicStructure> dataToBeFit;

        if (_externalPotential.has_value()) {
            JGAP_LOG_INFO("Subtracting external contributions");
            dataToBeFit = subtractExternalContribution(trainingData, _externalPotential.value());
        } else {
            dataToBeFit = std::vector(trainingData);
        }

        std::map<std::string, std::shared_ptr<Potential>> resultingPotentials;
        for (const auto &label: _fitOrder) {

            JGAP_LOG_INFO("Doing \"{}\" potential fit", label);
            resultingPotentials[label] = _fits[label] -> fit(dataToBeFit);
            JGAP_LOG_DEBUG(
                "Fitting finished for {}, resulting in : {}",
                label, resultingPotentials[label]->serialize().dump()
                );

            if (label != _fitOrder.back()) {
                dataToBeFit = subtractExternalContribution(dataToBeFit, resultingPotentials[label]);
            }
        }

        if (_externalPotential.has_value()) {
            // WARN: !! label reserved !!
            resultingPotentials["external"] = _externalPotential.value();
        }

        return std::make_shared<CompositePotential>(resultingPotentials);
    }

    std::vector<AtomicStructure> CompositePotentialFit::subtractExternalContribution(
        const std::vector<AtomicStructure> &originalData, const std::shared_ptr<Potential> &potential) {

        JGAP_LOG_INFO("Subtracting external potential contributions");

        std::vector dataToBeFit(originalData.begin(), originalData.end());
        NeighbourFinder::findNeighbours(dataToBeFit, potential->getCutoff().maxOverall());

        tbb::parallel_for_each(dataToBeFit.begin(), dataToBeFit.end(), [&](AtomicStructure& structure) -> void {
            structure.adjust(potential->predict(structure), true, false);
        });

        return dataToBeFit;
    }
}
