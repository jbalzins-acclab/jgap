#ifdef JGAP_TABULATIONDATA_HPP
#define JGAP_TABULATIONDATA_HPP

#include <map>
#include <ranges>
#include <set>

#include "io/log/CurrentLogger.hpp"

namespace jgap {

    struct TwoBodyTable {

    };

    class EamTabulationData {
    public:
        std::map<Species, Grid1d> densityGrids;
        std::map<ContributorReceiverSpecies, Grid1d> eamPairFunctionGrids;

        EamTabulationData(const std::optional<TabulationParams> &params) : _params(params) {}

        Grid1d& getOrMakeEnergyGrid(const Species& species) {
            if (!densityGrids.contains(species)) {
                if (!_params.has_value() || !_params->maxDensity.has_value() || !_params->minDensity.has_value()) {
                    JGAP_LOG_AND_THROW("Invalid EAM tabulation params -> can't make new density grids");
                }

                densityGrids[species] = Grid1d(
                    _params->nEamDensities,
                    (_params->maxDensity.value() -_params->minDensity.value())
                                / static_cast<double>(_params->nEamDensities),
                    _params->minDensity.value()
                    );
            }
            return densityGrids[species];
        }

        Grid1d getEnergyGridOrNull(const Species& species) const {
            if (!densityGrids.contains(species)) {
                return Grid1d(
                    _params->nEamDensities,
                    (_params->maxDensity.value() -_params->minDensity.value())
                                / static_cast<double>(_params->nEamDensities),
                    _params->minDensity.value()
                    );
            }
            return densityGrids.at(species);
        }

        Grid1d& getOrMakePairFunctionGrid(const ContributorReceiverSpecies& species) {
            if (!eamPairFunctionGrids.contains(species)) {
                if (!_params.has_value() || !_params->cutoff2b.has_value()) {
                    JGAP_LOG_AND_THROW("Invalid EAM tabulation params -> can't make std::pair function");
                }

                eamPairFunctionGrids[species] = Grid1d(
                    _params->n2b,
                    _params->cutoff2b.value() / static_cast<double>(_params->n2b),
                    0.0
                    );
            }
            return eamPairFunctionGrids[species];
        }

        Grid1d getPairFunctionGridOrNull(const ContributorReceiverSpecies& species) const {
            if (!eamPairFunctionGrids.contains(species)) {
                return Grid1d(
                    _params->n2b,
                    _params->cutoff2b.value() / static_cast<double>(_params->n2b),
                    0.0
                    );
            }
            return eamPairFunctionGrids.at(species);
        }

    private:
        const std::optional<TabulationParams>& _params;
    };

    class TabulationData {
    public:
        std::map<Species, Real> isolatedEnergies{};
        std::map<SpeciesPair, Grid1d> pairGrids{};
        std::vector<EamTabulationData> eamTabulationData{};
        std::map<SpeciesTriplet, Grid3d> tripletGrids{};

        TabulationData() : _params({}) {}
        TabulationData(const TabulationParams &params) : _params(params) {}

        Grid1d& getOrMake2bGrid(const SpeciesPair& species) {
            if (!pairGrids.contains(species)) {
                if (!_params.has_value() || !_params->cutoff2b.has_value()) {
                    JGAP_LOG_AND_THROW("Invalid 2b tabulation params -> can't 2b Grid");
                }

                pairGrids[species] = Grid1d(
                    _params->n2b,
                    _params->cutoff2b.value() / static_cast<double>(_params->n2b - 1),
                    0.0
                    );
            }
            return pairGrids[species];
        }

        Grid1d getOrNull(const SpeciesPair& species) const {
            if (!pairGrids.contains(species)) {
                return Grid1d(
                    _params->n2b,
                    _params->cutoff2b.value() / static_cast<double>(_params->n2b - 1),
                    0.0
                    );
            }
            return pairGrids.at(species);
        }

        Grid3d& getOrMake3bGrid(const SpeciesTriplet& species) {
            if (!tripletGrids.contains(species)) {
                if (!_params.has_value() || !_params->cutoff3b.has_value()) {
                    JGAP_LOG_AND_THROW("Invalid 3b tabulation params -> can't make 3b grid");
                }
                tripletGrids[species] = Grid3d(
                    _params->rMin3b.value(),
                    _params->n3bR,
                    (_params->cutoff3b.value() - _params->rMin3b.value()) / static_cast<double>(_params->n3bR - 1),
                    -1.0,
                    _params->n3bCosTheta,
                    2.0 / static_cast<double>(_params->n3bCosTheta - 1)
                    );
            }
            return tripletGrids[species];
        }

        EamTabulationData& newEamGrids() {
            eamTabulationData.emplace_back(_params);
            return eamTabulationData.back();
        }


    private:
        std::optional<TabulationParams> _params;
    };
}

#endif
