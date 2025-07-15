//
// Created by Jegors Balzins on 15.7.2025.
//

#include "core/fit/Tabulate.hpp"

#include <highfive/H5File.hpp>
#include <fstream>
#include <string>

#include "utils/AtomicNumbers.hpp"

using namespace std;

namespace jgap {

    void Tabulate::tabulate(const shared_ptr<Potential> &potential,
                            const nlohmann::json& params,
                            const string &outputFileNamePrefix) {

        const TabulationParams tabulationParams = parse(params);

        const TabulationData tabulationData = potential->tabulate(tabulationParams);

        if (tabulationData.eamTabulationData.size() > 1) {
            throw runtime_error("Multi EAM grid is not supported yet");
        }

        writeH5(potential, params, tabulationParams, tabulationData, outputFileNamePrefix);
        for (size_t index = 0; index < tabulationData.eamTabulationData.size(); index++) {
            writeEamFs(potential, params, tabulationParams,
                       tabulationData.eamTabulationData[index], outputFileNamePrefix,
                       tabulationData.pairEnergies, index == 0);
        }
    }

    void Tabulate::writeH5(const shared_ptr<Potential> &potential,
                           const nlohmann::json& params,
                           const TabulationParams &tabulationParams,
                           const TabulationData &tabulationData,
                           const string &outputFileNamePrefix) {

        HighFive::File tabgapFile(outputFileNamePrefix + ".tabgap.h5", HighFive::File::Overwrite);

        string comment1 = format("Tabulated jGAP: {}", potential->serialize().dump());
        tabgapFile.createDataSet<string>("comment1", HighFive::DataSpace::From(comment1)).write(comment1);
        string comment2 = format("Tabulation params: {}", params.dump());
        tabgapFile.createDataSet<string>("comment2", HighFive::DataSpace::From(comment2)).write(comment2);

        auto e0Group = tabgapFile.createGroup("e0");
        e0Group.createAttribute("Nelements", tabulationData.isolatedEnergies.size())
               .write(tabulationData.isolatedEnergies.size());
        for (const auto& [element, value] : tabulationData.isolatedEnergies) {
            e0Group.createAttribute(element, value).write(value);
        }

        if (tabulationData.eamTabulationData.empty()) {
            tabgapFile.createDataSet(
                "npots", vector{tabulationData.pairEnergies.size(), tabulationData.tripletEnergies.size()}
                );

            for (const auto& [speciesPair, energies] : tabulationData.pairEnergies) {
                auto pairGroup = tabgapFile.createGroup(format("{}-{}", speciesPair.first(), speciesPair.second()));

                pairGroup.createAttribute("element_i", speciesPair.first()).write(speciesPair.first());
                pairGroup.createAttribute("element_j", speciesPair.second()).write(speciesPair.second());

                pairGroup.createDataSet(
                    "grid_limits", vector{tabulationParams.grid2b[0], tabulationParams.grid2b.back()}
                    );

                pairGroup.createAttribute("N", tabulationParams.grid2b.size()).write(tabulationParams.grid2b.size());

                pairGroup.createDataSet("energies", energies);
            }

        } else {
            tabgapFile.createDataSet("npots", vector{0, tabulationData.tripletEnergies.size()});
        }

        for (const auto& [speciesTriplet, energies] : tabulationData.tripletEnergies) {
            auto tripletGroup = tabgapFile.createGroup(
                format("{}-{}-{}", speciesTriplet.root, speciesTriplet.nodes.first(), speciesTriplet.nodes.second())
                );
            tripletGroup.createAttribute("element_i", speciesTriplet.root).write(speciesTriplet.root);
            tripletGroup.createAttribute("element_j", speciesTriplet.nodes.first())
                        .write(speciesTriplet.nodes.first());
            tripletGroup.createAttribute("element_k", speciesTriplet.nodes.second())
                        .write(speciesTriplet.nodes.second());

            tripletGroup.createDataSet("grid_limits", vector{
                // WARN: convention sensitive - lowest in all first : highest in all last
                tabulationParams.grid3b[0][0], tabulationParams.grid3b[0][1], tabulationParams.grid3b[0][2],
                tabulationParams.grid3b.back()[0], tabulationParams.grid3b.back()[1], tabulationParams.grid3b.back()[2],
            });

            tripletGroup.createDataSet("N", vector{
                params.value("n3b_r", 80),
                params.value("n3b_r", 80),
                params.value("n3b_angle", 80)
            });

            tripletGroup.createDataSet("energies", energies);
        }

        tabgapFile.flush();
    }

    // TODO: looks monstrous
    void Tabulate::writeEamFs(const shared_ptr<Potential> &potential,
                              const nlohmann::json& params,
                              const TabulationParams &tabulationParams,
                              const EamTabulationData& data,
                              const string &outputFileNamePrefix,
                              const map<SpeciesPair, vector<double>> &pairEnergies,
                              const bool writePairEnergies) {

        ofstream eamFsFile(outputFileNamePrefix + ".eam.fs");

        // Lines 1–3: Comments/metadata.
        eamFsFile << "# UNITS: metal" << endl;
        eamFsFile << "# Tabulated jGAP: <<" << potential->serialize().dump()
                  << ">> with params: <<" << params.dump() << endl;
        eamFsFile << "# pair_style eam/fs" << endl;

        // Line 4: Number of elements (N) followed by each element’s symbol
        eamFsFile << data.embeddingEnergies.size() << " ";
        vector<Species> elements; // keep indexing consistent
        for (const auto& species : data.embeddingEnergies | views::keys) {
            elements.emplace_back(species);
            eamFsFile << elements.back() << " ";
        }
        eamFsFile << endl;

        // Line 5: Nrho, drho, Nr, dr, cutoff
        eamFsFile << tabulationParams.nDensities << " ";
        eamFsFile << data.rhoMax / static_cast<double>(tabulationParams.nDensities) << " ";
        eamFsFile << tabulationParams.grid2b.size() << " ";
        eamFsFile << tabulationParams.grid2b[1] << " ";
        eamFsFile << tabulationParams.grid2b.back() << endl;

        // Per-element Sections:
        /*
            * Line: atomic number, mass, lattice constant, lattice type (e.g., fcc, bcc)
            * Embedding function F_\beta(\rho): Nrho tabulated values
            * Density functions \rho_{\alpha\beta}(r): For each element α (total N curves, each with Nr points)
         */
        for (const Species& species1 : elements) {
            eamFsFile << Z_default[species1] << " ";
            eamFsFile << mass_default[species1] << " 1.0 ZZZ" << endl;

            for (const double& energy : data.embeddingEnergies.at(species1)) {
                eamFsFile << energy << endl;
            }

            for (const Species& species2 : elements) {
                for (const double& density : data.eamDensities.at(OrderedSpeciesPair{species1, species2})) {
                    eamFsFile << density << endl;
                }
            }
        }

        /*
        *Pair Potential Tables (for all i ≥ j pairs):
            *Tabulated \phi_{\alpha\beta}(r) values for each unique pair (symmetry is exploited),
            *listing Nr points per interaction
         */

        for (size_t i = 0; i < elements.size(); i++) {
            for (size_t j = i; j < elements.size(); j++) {
                for (const double& energy: pairEnergies.at({elements[i], elements[j]})) {
                    eamFsFile << (writePairEnergies ? energy : 0) << endl;
                }
            }
        }

        eamFsFile.flush();
        eamFsFile.close();
    }

    TabulationParams Tabulate::parse(const nlohmann::json &params) {

        TabulationParams result{};
        result.grid2b = make2bGrid(params);
        result.grid3b = make3bGrid(params);
        result.nDensities = params.value("nDensities", 1000);

        result.species = {};
        for (const auto& element: params["species"]) {
            result.species.emplace_back(element.get<string>());
        }

        return result;
    }

    vector<double> Tabulate::make2bGrid(const nlohmann::json &params) {

        const size_t n2b = params.value("n2b", 1000);
        // TODO: add option to customize?
        const double rMin2b = 0.0; //params.value("r_min_2b", 0.0);
        const double rMax2b = params["r_max_2b"];
        const double dr_2b = (rMax2b - rMin2b) / static_cast<double>(n2b - 1);

        vector<double> grid2b{};
        for (size_t i = 0; i < n2b; i++) {
            grid2b.push_back(rMin2b + static_cast<double>(i) * dr_2b);
        }
        grid2b[n2b - 1] = rMax2b; // just in case
        return grid2b;
    }

    vector<Vector3> Tabulate::make3bGrid(const nlohmann::json &params) {

        const double rMin3b = params.value("r_min_3b", 0.1);
        const double rMax3b = params["r_max_3b"];
        const size_t n3b_r = params.value("n3b_r", 80);
        const double dr_3b = (rMax3b - rMin3b) / static_cast<double>(n3b_r - 1);

        const size_t n3b_angle = params.value("n3b_angle", 80);
        const double angleStep = 2.0/*cos: from -1 to 1*/ / static_cast<double>(n3b_angle - 1);

        vector<Vector3> result{};
        for (size_t i = 0; i < n3b_r; i++) {
            double r1 = rMin3b + static_cast<double>(i) * dr_3b;

            for (size_t j = 0; j < n3b_r; j++) {
                double r2 = rMin3b + static_cast<double>(j) * dr_3b;

                for (size_t k = 0; k < n3b_angle; k++) {
                    double angle = angleStep * static_cast<double>(k) + angleStep;
                    result.push_back(Vector3{r1, r2, angle});
                }
            }
        }

        return result;
    }


}
