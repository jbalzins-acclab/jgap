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

        CurrentLogger::get()->debug("Starting tabulation");
        const TabulationData tabulationData = potential->tabulate(tabulationParams);
        CurrentLogger::get()->debug("Finished tabulation");

        CurrentLogger::get()->debug("Saving H5");
        writeH5(potential, params, tabulationParams, tabulationData, outputFileNamePrefix);

        for (size_t index = 0; index < tabulationData.eamTabulationData.size(); index++) {
            CurrentLogger::get()->debug(format("Saving eam.fs #{}", index));
            writeEamFs(potential, params, tabulationParams,
                       tabulationData.eamTabulationData[index], outputFileNamePrefix,
                       tabulationData.pairEnergies, index);
        }
    }

    void Tabulate::writeH5(const shared_ptr<Potential> &potential,
                           const nlohmann::json& params,
                           const TabulationParams &tabulationParams,
                           const TabulationData &tabulationData,
                           const string &outputFileNamePrefix) {

        HighFive::File tabGapFile(outputFileNamePrefix + ".tabgap.h5", HighFive::File::Overwrite);

        const string comment1 = format("Tabulated jGAP: {}", potential->serialize().dump());
        tabGapFile.createDataSet<string>("comment1", HighFive::DataSpace::From(comment1)).write(comment1);
        const string comment2 = format("Tabulation params: {}", params.dump());
        tabGapFile.createDataSet<string>("comment2", HighFive::DataSpace::From(comment2)).write(comment2);

        auto e0Group = tabGapFile.createGroup("e0");
        e0Group.createAttribute("Nelements", tabulationData.isolatedEnergies.size())
               .write(tabulationData.isolatedEnergies.size());
        for (const auto& [element, value]: tabulationData.isolatedEnergies) {
            e0Group.createAttribute(element, value).write(value);
        }

        if (tabulationData.eamTabulationData.empty()) {
            tabGapFile.createDataSet(
                "npots", vector{tabulationData.pairEnergies.size(), tabulationData.tripletEnergies.size()}
                );

            for (const auto& [speciesPair, energies] : tabulationData.pairEnergies) {
                auto pairGroup = tabGapFile.createGroup(format("{}-{}", speciesPair.first(), speciesPair.second()));

                pairGroup.createAttribute("element_i", speciesPair.first()).write(speciesPair.first());
                pairGroup.createAttribute("element_j", speciesPair.second()).write(speciesPair.second());

                pairGroup.createDataSet(
                    "grid_limits", vector{tabulationParams.grid2b[0], tabulationParams.grid2b.back()}
                    );

                pairGroup.createAttribute("N", tabulationParams.grid2b.size()).write(tabulationParams.grid2b.size());

                auto splineCoefficients = toSplineCoefficients(
                    energies, tabulationParams.grid2b[1] - tabulationParams.grid2b[0]
                    );
                pairGroup.createDataSet("energies", splineCoefficients);
            }

        } else {
            tabGapFile.createDataSet("npots", vector{0, tabulationData.tripletEnergies.size()});
        }

        for (const auto& [speciesTriplet, energies] : tabulationData.tripletEnergies) {
            auto tripletGroup = tabGapFile.createGroup(
                format("{}-{}-{}", speciesTriplet.root, speciesTriplet.nodes.first(), speciesTriplet.nodes.second())
                );
            tripletGroup.createAttribute("element_i", speciesTriplet.root).write(speciesTriplet.root);
            tripletGroup.createAttribute("element_j", speciesTriplet.nodes.first())
                        .write(speciesTriplet.nodes.first());
            tripletGroup.createAttribute("element_k", speciesTriplet.nodes.second())
                        .write(speciesTriplet.nodes.second());

            tripletGroup.createDataSet("grid_limits", vector{
                // WARN: convention sensitive - lowest in all first : highest in all last
                tabulationParams.grid3b[0][0][0].x,
                tabulationParams.grid3b[0][0][0].y,
                tabulationParams.grid3b[0][0][0].z,
                tabulationParams.grid3b.back().back().back().x, // :)
                tabulationParams.grid3b.back().back().back().y,
                tabulationParams.grid3b.back().back().back().z,
            });

            tripletGroup.createDataSet("N", vector{
                params.value("n3b_r", 80),
                params.value("n3b_r", 80),
                params.value("n3b_angle", 80)
            });

            auto splineCoefficients = toSplineCoefficients(
                energies, tabulationParams.grid3b[1][1][1] - tabulationParams.grid3b[0][0][0]
                );
            tripletGroup.createDataSet("energies", splineCoefficients);
        }

        tabGapFile.flush();
    }

    // TODO: looks monstrous
    void Tabulate::writeEamFs(const shared_ptr<Potential> &potential,
                              const nlohmann::json& params,
                              const TabulationParams &tabulationParams,
                              const EamTabulationData& data,
                              const string &outputFileNamePrefix,
                              const map<SpeciesPair, vector<double>> &pairEnergies,
                              const size_t index) {

        const string filename = outputFileNamePrefix + (index != 0 ? "#" + to_string(index) : "") + ".eam.fs";

        ofstream eamFsFile(filename);
        if (!eamFsFile.is_open()) {
            CurrentLogger::get()->error("Could not open " + filename, true);
        }
        eamFsFile << fixed << setprecision(17);

        // Lines 1–3: Comments/metadata.
        eamFsFile << "# UNITS: metal" << endl;
        eamFsFile << "# Tabulated jGAP"  << endl;
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
        eamFsFile << data.maxDensity / static_cast<double>(tabulationParams.nDensities - 1) << " ";
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
            eamFsFile << static_cast<size_t>(Z_default[species1]) << " ";
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
            for (size_t j = 0; j < elements.size(); j++) {
                if (i < j) {
                    continue; // AAAAAAAA why not i <= j????????
                }
                //cout << i << " " << j << endl;
                //cout << elements[i] << " " << elements[j] << endl;
                const auto& energies = pairEnergies.at({elements[i], elements[j]});
                for (size_t k = 0; k < energies.size(); k++) {
                    eamFsFile << (index == 0 ? energies[k]*tabulationParams.grid2b[k] : 0) << endl;
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
        result.nDensities = params.value("n_densities", 1000);
        if (params.contains("max_eam_density")) {
            result.maxDensity = params["max_eam_density"];
        }

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
        grid2b.reserve(n2b);
        for (size_t i = 0; i < n2b; i++) {
            grid2b.push_back(rMin2b + static_cast<double>(i) * dr_2b);
        }
        grid2b[n2b - 1] = rMax2b; // just in case
        return grid2b;
    }

    vector<vector<vector<Vector3>>> Tabulate::make3bGrid(const nlohmann::json &params) {

        const double rMin3b = params.value("r_min_3b", 0.1);
        const double rMax3b = params["r_max_3b"];
        const size_t n3b_r = params.value("n3b_r", 80);
        const double dr_3b = (rMax3b - rMin3b) / static_cast<double>(n3b_r - 1);

        const size_t n3b_angle = params.value("n3b_angle", 80);
        const double angleStep = 2.0/*cos: from -1 to 1*/ / static_cast<double>(n3b_angle - 1);

        vector result(n3b_r, vector(n3b_r, vector(n3b_angle, Vector3(0, 0, 0))));
        for (size_t i = 0; i < n3b_r; i++) {
            double r1 = rMin3b + static_cast<double>(i) * dr_3b;

            for (size_t j = 0; j < n3b_r; j++) {
                double r2 = rMin3b + static_cast<double>(j) * dr_3b;

                for (size_t k = 0; k < n3b_angle; k++) {
                    double angle = -1.0 + angleStep * static_cast<double>(k);
                    result[i][j][k] = {r1, r2, angle};
                }
            }
        }

        return result;
    }

    vector<double> Tabulate::toSplineCoefficients(const vector<double>& energies, const double spacing) {
        static constexpr array basis = {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0};

        const size_t nCoefficients = energies.size() + 2;

        const double inverseSpacing = 1.0 / spacing;
        const double inverseSpacingSq = inverseSpacing * inverseSpacing;
        vector<array<double, 4>> bands(nCoefficients);
        bands[0] = {inverseSpacingSq, -2.0 * inverseSpacingSq, inverseSpacingSq, 0.0};
        //bands[0] = {1.0, -2.0, 1.0, 0.0};
        bands[nCoefficients-1] = {inverseSpacingSq, -2.0 * inverseSpacingSq, inverseSpacingSq, 0.0};

        for (size_t i = 1; i < nCoefficients-1; i++) {
            bands[i][0] = basis[0];
            bands[i][1] = basis[1];
            bands[i][2] = basis[2];
            bands[i][3] = energies[i - 1];
        }

        bands[0][1] /= bands[0][0];
        bands[0][2] /= bands[0][0];
        bands[0][3] /= bands[0][0];
        bands[0][0] = 1.0;
        bands[1][1] -= bands[1][0] * bands[0][1];
        bands[1][2] -= bands[1][0] * bands[0][2];
        bands[1][3] -= bands[1][0] * bands[0][3];
        bands[0][0] = 0;
        bands[1][2] /= bands[1][1];
        bands[1][3] /= bands[1][1];
        bands[1][1] = 1.0;

        for (size_t i = 2; i < nCoefficients-1; i++) {
            bands[i][1] -= bands[i][0] * bands[i-1][2];
            bands[i][3] -= bands[i][0] * bands[i-1][3];
            bands[i][2] /= bands[i][1];
            bands[i][3] /= bands[i][1];
            bands[i][0] = 0.0;
            bands[i][1] = 1.0;
        }

        bands[nCoefficients-1][1] -= bands[nCoefficients-1][0] * bands[nCoefficients-3][2];
        bands[nCoefficients-1][3] -= bands[nCoefficients-1][0] * bands[nCoefficients-3][3];
        bands[nCoefficients-1][2] -= bands[nCoefficients-1][1] * bands[nCoefficients-2][2];
        bands[nCoefficients-1][3] -= bands[nCoefficients-1][1] * bands[nCoefficients-2][3];
        bands[nCoefficients-1][3] /= bands[nCoefficients-1][2];
        bands[nCoefficients-1][2] = 1.0;

        vector coefficients(nCoefficients, 0.0);
        coefficients[nCoefficients-1] = bands[nCoefficients-1][3];
        for (size_t i = nCoefficients-2; i > 0; i--) {
            coefficients[i] = bands[i][3] - bands[i][2] * coefficients[i+1];
        }
        coefficients[0] = bands[0][3] - bands[0][1] * coefficients[1] - bands[0][2] * coefficients[2];

        return coefficients;
    }

    vector<double> Tabulate::toSplineCoefficients(const vector<vector<vector<double>>>& energies,
                                                   const Vector3 &spacing) {

        const array nCoefficients{energies.size()+2, energies[0].size()+2, energies[0][0].size()+2};
        vector coefficients(nCoefficients[0], vector(nCoefficients[1], vector(nCoefficients[2], 0.0)));

        // For convenience starting with angle
        for (size_t i = 0; i < energies.size(); i++) {
            for (size_t j = 0; j < energies[i].size(); j++) {
                for (size_t k = 0; k < energies[i][j].size(); k++) {
                    coefficients[i+1][j+1] = toSplineCoefficients(energies[i][j], spacing.z);
                }
            }
        }

        for (size_t i = 0; i < energies.size(); i++) {
            for (size_t k = 0; k < nCoefficients[2]; k++) {
                vector coeffs_i_k(energies[i].size(), 0.0);
                for (size_t j = 0; j < energies[i].size(); j++) {
                    coeffs_i_k[j] = coefficients[i+1][j+1][k];
                }
                vector coeffs_i_k_full = toSplineCoefficients(coeffs_i_k, spacing.y);
                for (size_t j = 0; j < nCoefficients[1]; j++) {
                    coefficients[i+1][j][k] = coeffs_i_k_full[j];
                }
            }
        }

        for (size_t j = 0; j < nCoefficients[1]; j++) {
            for (size_t k = 0; k < nCoefficients[2]; k++) {
                vector coeffs_j_k(energies.size(), 0.0);
                for (size_t i = 0; i < energies.size(); i++) {
                    coeffs_j_k[i] = coefficients[i+1][j][k];
                }
                vector coeffs_j_k_full = toSplineCoefficients(coeffs_j_k, spacing.x);
                for (size_t i = 0; i < nCoefficients[0]; i++) {
                    coefficients[i][j][k] = coeffs_j_k_full[i];
                }
            }
        }

        vector<double> coefficientsFlat;
        coefficientsFlat.reserve(nCoefficients[0] * nCoefficients[1] * nCoefficients[2]);
        for (size_t i = 0; i < nCoefficients[0]; i++) {
            for (size_t j = 0; j < nCoefficients[1]; j++) {
                for (size_t k = 0; k < nCoefficients[2]; k++) {
                    coefficientsFlat.push_back(coefficients[i][j][k]);
                }
            }
        }

        return coefficientsFlat;
    }
}
