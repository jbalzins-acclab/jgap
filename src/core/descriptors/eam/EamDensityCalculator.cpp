#include "core/descriptors/eam/EamDensityCalculator.hpp"

#include "core/descriptors/eam/pair_functions/CoscutoffPairFunction.hpp"
#include "../../../../include/core/descriptors/eam/pair_functions/species/EamSpeciesPairFunction.hpp"
#include "core/descriptors/eam/pair_functions/FSGenPairFunction.hpp"
#include "../../../../include/core/descriptors/eam/pair_functions/species/FsgenSpeciesPairFunction.hpp"
#include "../../../../include/core/descriptors/eam/pair_functions/species/FssymSpeciesPairFunction.hpp"
#include "core/descriptors/eam/pair_functions/PolycutoffPairFunction.hpp"

namespace jgap {

    EamDensityCalculator::EamDensityCalculator(EamDescriptorParams params) {

        _cutoff = 0;

        _defaultPairFunction = nullptr;
        if (params.defaultDensityCalculationParams.has_value()) {
            _cutoff = params.defaultDensityCalculationParams.value().cutoff;
            _defaultPairFunction = selectPairFunction(params.defaultDensityCalculationParams.value(),{});
        } else {
            if (params.densityCalculationParamsPerSpecies.empty()) {
                throw runtime_error("EAM density calculation not parametrized!");
            }
        }

        _pairFunctions = {};
        for (auto &[orderedSpeciesPair, paramsPerSpecies]: params.densityCalculationParamsPerSpecies) {
            if (paramsPerSpecies.cutoff > _cutoff) _cutoff = paramsPerSpecies.cutoff;
            _pairFunctions[orderedSpeciesPair] = selectPairFunction(paramsPerSpecies, orderedSpeciesPair);
        }
    }

    EamKernelIndex EamDensityCalculator::calculate(const AtomicStructure &structure) const {

        EamKernelIndex result{};

        for (size_t atomIdx = 0; atomIdx < structure.atoms.size(); atomIdx++) {

            double totalDensity = 0;
            vector<pair<NeighbourData, double>> densityDerivatives;

            for (NeighbourData neighbour : structure.atoms[atomIdx].neighbours.value()) {
                if (neighbour.distance > _cutoff) continue;

                pair orderedSpeciesPair = {
                    structure.atoms[neighbour.index].species, structure.atoms[atomIdx].species
                };

                if (!_pairFunctions.contains(orderedSpeciesPair)) {
                    if (_defaultPairFunction == nullptr) continue;

                    totalDensity += _defaultPairFunction->evaluate(neighbour.distance);
                    densityDerivatives.emplace_back(
                        neighbour,
                        _defaultPairFunction->differentiate(neighbour.distance)
                    );
                }
            }

            result.push_back({totalDensity, densityDerivatives});
        }

        return result;
    }

    shared_ptr<EamPairFunction> EamDensityCalculator::selectPairFunction(EamDensityCalculationParams params,
                                                                         optional<OrderedSpeciesPair> speciesPair) {
        shared_ptr<EamPairFunction> result;

        switch (params.pairFunction) {
            case EamDensityCalculationParams::PairFunction::FSGEN:
                result = make_shared<FSGenPairFunction>(params.cutoff, params.degree.value());
                break;
            case EamDensityCalculationParams::PairFunction::COSCUTOFF:
                result = make_shared<CoscutoffPairFunction>(params.cutoff, params.rmin.value());
                break;
            case EamDensityCalculationParams::PairFunction::POLYCUTOFF:
                result = make_shared<PolycutoffPairFunction>(params.cutoff, params.rmin.value());
                break;
            default:
                throw runtime_error("Pair function not specified");
        }

        if (speciesPair.has_value()) {
            switch (params.speciesUse) {
                case EamDensityCalculationParams::SpeciesUse::EAM:
                    result = make_shared<EamSpeciesPairFunction>(result, speciesPair->second);
                    break;
                case EamDensityCalculationParams::SpeciesUse::FSGEN:
                    result = make_shared<FsgenSpeciesPairFunction>(result, speciesPair->first, speciesPair->second);
                    break;
                case EamDensityCalculationParams::SpeciesUse::FSSYM:
                    result = make_shared<FssymSpeciesPairFunction>(result, speciesPair->first, speciesPair->second);
                    break;
                case EamDensityCalculationParams::SpeciesUse::BLIND:
                    // e.g. all but one pair is blind
                    break;
                default:
                    throw runtime_error("Species use in EAM not specified correcty for species pair - "
                                        + speciesPair->first + speciesPair->second);
            }
        } else {

        }

        return result;
    }
}
