#include "core/fit/CompositeFit.hpp"
#include "core/potentials/CompositePotential.hpp"

namespace jgap {
    CompositePotentialFit::CompositePotentialFit(const nlohmann::json &params) {
        if (params.contains("external")) {
            _externalPotential = ParserRegistry<Potential>::get(params["external"]);
        } else {
            _externalPotential = {};
        }

        _fits = {};
        for (const auto &[label, fitParams]: params["fits"].items()) {
            CurrentLogger::get()->info("Picking fitting logic for " + label);
            _fits[label] = ParserRegistry<Fit>::get(fitParams);
        }

        if (params.contains("fit_order")) {
            for (const auto& label: params["fit_order"]) {
                _fitOrder.push_back(label);
            }
            if (_fitOrder.size() != _fits.size()) {
                CurrentLogger::get()->warn("Fit order size mismatches number of fit params");
            }
        } else {
            if (_fits.size() == 1) {
                _fitOrder = { _fits.begin()->first };
            } else if (_fits.size() == 2 && _fits.contains("isolated_atom")) {
                _fitOrder = vector(views::keys(_fits).begin(), views::keys(_fits).end());
                if (_fitOrder[0] != "isolated_atom") {
                    swap(_fitOrder[0], _fitOrder[1]);
                }
            } else {
                CurrentLogger::get()->error("Fitting order not specified", true);
            }
        }
    }

    shared_ptr<Potential> CompositePotentialFit::fit(const vector<AtomicStructure> &trainingData) {
        vector<AtomicStructure> dataToBeFit;

        if (_externalPotential.has_value()) {
            CurrentLogger::get()->info("Subtracting external contributions");
            dataToBeFit = subtractExternalContribution(trainingData, _externalPotential.value());
        } else {
            dataToBeFit = vector(trainingData);
        }

        map<string, shared_ptr<Potential>> resultingPotentials;
        for (const auto &label: _fitOrder) {

            CurrentLogger::get()->info(format("Doing \"{}\" potential fit", label));
            resultingPotentials[label] = _fits[label] -> fit(dataToBeFit);
            CurrentLogger::get()->debug(format(
                "Fitting finished for {}, resulting in : {}",
                label, resultingPotentials[label]->serialize().dump()
                ));

            if (label != _fitOrder.back()) {
                dataToBeFit = subtractExternalContribution(dataToBeFit, resultingPotentials[label]);
            }
        }

        /*
        nlohmann::json potentialsSerialized{};
        if (_externalPotential.has_value()) {
            potentialsSerialized["external"] = _externalPotential.value()->serialize();
            potentialsSerialized["external"]["type"] = _externalPotential.value()->getType();
        }

        for (const auto &[label, fittedPotential]: fittedPotentials) {
            potentialsSerialized[label] = fittedPotential->serialize();
            potentialsSerialized[label]["type"] = fittedPotential->getType();
        }
        */

        if (_externalPotential.has_value()) {
            // WARN: !! label reserved !!
            resultingPotentials["external"] = _externalPotential.value();
        }

        return make_shared<CompositePotential>(resultingPotentials);
    }

    vector<AtomicStructure> CompositePotentialFit::subtractExternalContribution(
        const vector<AtomicStructure> &originalData, const shared_ptr<Potential> &potential) {

        CurrentLogger::get()->info("Subtracting external potential contributions");

        vector dataToBeFit(originalData.begin(), originalData.end());
        NeighbourFinder::findNeighbours(dataToBeFit, potential->getCutoff());

        for (auto& structure: dataToBeFit) {
            structure.adjust(potential->predict(structure), true, false);
        }

        return dataToBeFit;
    }
}
