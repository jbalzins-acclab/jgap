#ifndef JGAP_TABULATIONDATA_HPP
#define JGAP_TABULATIONDATA_HPP

#include <map>

#include "Grid1d.hpp"
#include "Grid3d.hpp"
#include "SpeciesData.hpp"
#include "data/Vector3.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {

    struct TabulationParams {
        size_t n2b;
        double cutoff2b;

        size_t nEamDensities;
        double minDensity, maxDensity;

        size_t n3bR, n3bCosTheta;
        double rMin3b, cutoff3b;
    };

    class EamData {
    public:
        const TabulationParams& params;

        map<Species, Grid1d> embeddingEnergies;
        map<ContributorReceiverSpecies, Grid1d> eamPairFunctions;

        EamData(const TabulationParams &params) : params(params) {}

        Grid1d& getEnergyGrid(const Species& species) {
            if (!embeddingEnergies.contains(species)) {
                embeddingEnergies[species] = Grid1d(
                    params.nEamDensities,
                    (params.maxDensity - params.minDensity) / static_cast<double>(params.nEamDensities),
                    params.minDensity
                    );
            }
            return embeddingEnergies[species];
        }

        Grid1d& getPairFunctionGrid(const ContributorReceiverSpecies& species) {
            if (!eamPairFunctions.contains(species)) {
                eamPairFunctions[species] = Grid1d(
                    params.n2b,
                    params.cutoff2b / static_cast<double>(params.n2b),
                    0.0
                    );
            }
            return eamPairFunctions[species];
        }
    };

    class TabulationData {
    public:
        TabulationParams params;

        map<Species, double> isolatedEnergies{};
        map<SpeciesPair, Grid1d> pairEnergies{};
        vector<EamData> eamTabulationData{};
        map<SpeciesTriplet, Grid3d> tripletEnergies{};

        TabulationData(const TabulationParams &params) : params(params) {}

        Grid1d& get2bGrid(const SpeciesPair& species) {
            if (!pairEnergies.contains(species)) {
                pairEnergies[species] = Grid1d(
                    params.n2b,
                    params.cutoff2b / static_cast<double>(params.n2b - 1),
                    0.0
                    );
            }
            return pairEnergies[species];
        }

        Grid3d& get3bGrid(const SpeciesTriplet& species) {
            if (!tripletEnergies.contains(species)) {
                tripletEnergies[species] = Grid3d(
                    params.rMin3b,
                    params.n3bR,
                    (params.cutoff3b - params.rMin3b) / static_cast<double>(params.n3bR - 1),
                    0.0,
                    params.n3bCosTheta,
                    2.0 / static_cast<double>(params.n3bCosTheta - 1)
                    );
            }
            return tripletEnergies[species];
        }

        EamData& newEamGrids() {
            eamTabulationData.emplace_back(params);
            return eamTabulationData.back();
        }
    };
}

#endif
