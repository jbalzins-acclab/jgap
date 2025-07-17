//
// Created by Jegors Balzins on 15.7.2025.
//

#ifndef TABULATIONDATA_HPP
#define TABULATIONDATA_HPP

#include <map>
#include "data/BasicDataTypes.hpp"

namespace jgap {

    struct TabulationParams {
        // TODO: rho_max
        vector<Species> species{};
        size_t nDensities{};
        vector<double> grid2b{}; // nRho = grid2b.size() TODO ?
        vector<vector<vector<Vector3>>> grid3b{};
    };

    struct EamTabulationData {
        // rhoMin = 0 always!!
        double rhoMax;
        map<Species, vector<double>> embeddingEnergies;
        map<OrderedSpeciesPair, vector<double>> eamDensities; // per grid2b
    };

    struct TabulationData {

        map<Species, double> isolatedEnergies{};
        map<SpeciesPair, vector<double>> pairEnergies{};
        map<SpeciesTriplet, vector<vector<vector<double>>>> tripletEnergies{};

        vector<EamTabulationData> eamTabulationData{};

        TabulationData operator+(const TabulationData &other) const {
            TabulationData result = *this;

            for (const auto &[species, energy]: other.isolatedEnergies) {
                if (!result.isolatedEnergies.contains(species)) {
                    result.isolatedEnergies[species] = energy;
                } else {
                    result.isolatedEnergies[species] += energy;
                }
            }

            for (const auto &[speciesPair, energies]: other.pairEnergies) {
                if (!result.pairEnergies.contains(speciesPair)) {
                    result.pairEnergies[speciesPair] = energies;
                } else {
                    if (result.pairEnergies[speciesPair].size() != energies.size()) {
                        CurrentLogger::get()->error("Pair energy table size mismatch", true);
                    }
                    for (size_t i = 0; i < energies.size(); i++) {
                        result.pairEnergies[speciesPair][i] += energies[i];
                    }
                }
            }

            for (const auto &[speciesTriplet, energies]: other.tripletEnergies) {
                if (!result.tripletEnergies.contains(speciesTriplet)) {
                    result.tripletEnergies[speciesTriplet] = energies;
                } else {
                    if (result.tripletEnergies[speciesTriplet].size() != energies.size()) {
                        CurrentLogger::get()->error("Triplet energy table size mismatch", true);
                    }
                    for (size_t i = 0; i < energies.size(); i++) {
                        for (size_t j = 0; j < energies[i].size(); j++) {
                            for (size_t k = 0; k < energies[i][j].size(); k++) {
                                result.tripletEnergies[speciesTriplet][i][j][k] += energies[i][j][k];
                            }
                        }
                    }
                }
            }

            result.eamTabulationData = eamTabulationData;
            result.eamTabulationData.insert(result.eamTabulationData.end(),
                                            other.eamTabulationData.begin(),
                                            other.eamTabulationData.end());

            return result;
        }
    };
}

#endif //TABULATIONDATA_HPP
