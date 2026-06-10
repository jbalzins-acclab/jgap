#if false
#include "io/tabgap/TabGapIO.hpp"

#include <fstream>
#include <sstream>
#include <algorithm>
#include <highfive/H5File.hpp>


namespace jgap {

    TabulationData TabGapIO::read(const FileNames& fileNames) {

        TabulationData result;

        bool didReadH5 = false;

        for (const auto& fileName: fileNames) {
            JGAP_LOG_DEBUG("Reading file {}", fileName);
            if (fileName.ends_with(".h5")) {
                if (didReadH5) {
                    JGAP_LOG_AND_THROW(
                        "Multiple .h5 files detected when reading single tabulation data: {}",
                        vectorToString(fileNames)
                        );
                }
                readH5(fileName, result);
                didReadH5 = true;
            } else if (fileName.ends_with(".eam.fs")) {
                readEamFs(fileName, result);
            } else {
                JGAP_LOG_AND_THROW("Could not recognize file's extension: {}", fileName);
            }
        }

        return result;
    }

    FileNames TabGapIO::write(const TabulationData &valuesTables,
                                   const TabulationData &splineTables,
                                   std::optional<std::string> outputFileNamePrefix) {

        if (!outputFileNamePrefix.has_value()) {
            outputFileNamePrefix = generateFileNamePrefix(splineTables);
        }

        FileNames resultingFileNames;
        JGAP_LOG_DEBUG("Saving H5");
        resultingFileNames.push_back(writeH5(valuesTables, splineTables, outputFileNamePrefix.value()));

        for (size_t index = 0; index < valuesTables.eamTabulationData.size(); index++) {
            JGAP_LOG_DEBUG("Saving eam.fs #{}", index);
            writeEamFs(valuesTables, index, outputFileNamePrefix.value());
        }
        return resultingFileNames;
    }

    std::string TabGapIO::generateFileNamePrefix(const TabulationData &table) {
        std::string speciesStr;
        for (const std::string& species: table.allSpecies()) {
            speciesStr += species;
        }

        return speciesStr + uniqueStamp();
    }

    std::string TabGapIO::writeH5(const TabulationData &valuesTables,
                             const TabulationData &splineTables,
                             const std::string &outputFileNamePrefix) {

        // REMEMBER(!): Spline grid written, but with original grid specifications
        HighFive::File tabGapFile(outputFileNamePrefix + ".tabgap.h5", HighFive::File::Overwrite);

        const std::string comment1 = "UNITS: metal";
        tabGapFile.createDataSet<std::string>("comment1", HighFive::DataSpace::From(comment1)).write(comment1);
        const std::string comment2 = "pair_style tabgap";
        tabGapFile.createDataSet<std::string>("comment2", HighFive::DataSpace::From(comment2)).write(comment2);

        auto e0Group = tabGapFile.createGroup("e0");
        e0Group.createAttribute("Nelements", splineTables.isolatedEnergies.size())
               .write(splineTables.isolatedEnergies.size());
        for (const auto& [element, value]: splineTables.isolatedEnergies) {
            e0Group.createAttribute(element, value).write(value);
        }

        if (splineTables.eamTabulationData.empty()) {
            tabGapFile.createDataSet(
                "npots", std::vector{splineTables.pairGrids.size(), splineTables.tripletGrids.size()}
                );

            for (const auto& [speciesPair, coeffs2b] : splineTables.pairGrids) {
                auto pairGroup = tabGapFile.createGroup(format("{}-{}", speciesPair.first(), speciesPair.second()));

                pairGroup.createAttribute("element_i", speciesPair.first()).write(speciesPair.first());
                pairGroup.createAttribute("element_j", speciesPair.second()).write(speciesPair.second());

                const auto& pairEnergies = valuesTables.pairGrids.at(speciesPair);
                pairGroup.createDataSet( /*Original grid, not spline coeffs*/
                    "grid_limits", std::vector{pairEnergies.origin, pairEnergies.cutoff()}
                );

                pairGroup.createAttribute("N", pairEnergies.size())
                         .write(pairEnergies.size()/*Original grid, not spline coeffs*/);
                pairGroup.createDataSet("energies", coeffs2b.data);
            }

        } else {
            tabGapFile.createDataSet("npots", std::vector{0, splineTables.tripletGrids.size()});
        }

        for (const auto& [speciesTriplet, coeffs] : splineTables.tripletGrids) {
            auto tripletGroup = tabGapFile.createGroup(
                format("{}-{}-{}", speciesTriplet.root, speciesTriplet.nodes.first(), speciesTriplet.nodes.second())
                );
            tripletGroup.createAttribute("element_i", speciesTriplet.root).write(speciesTriplet.root);
            tripletGroup.createAttribute("element_j", speciesTriplet.nodes.first())
                        .write(speciesTriplet.nodes.first());
            tripletGroup.createAttribute("element_k", speciesTriplet.nodes.second())
                        .write(speciesTriplet.nodes.second());

            const auto& tripletEnergies = valuesTables.tripletGrids.at(speciesTriplet);

            tripletGroup.createDataSet("grid_limits", std::vector{
                tripletEnergies.origin.x, tripletEnergies.origin.y, -1.0,
                tripletEnergies.cutoff(), tripletEnergies.cutoff(), 1.0
            });

            tripletGroup.createDataSet("N", std::vector{
                tripletEnergies.nR, tripletEnergies.nR, tripletEnergies.nAngular
            });

            tripletGroup.createDataSet("energies", coeffs.dataFlat);
        }

        tabGapFile.flush();
        return tabGapFile.getName();
    }

    std::string TabGapIO::writeEamFs(const TabulationData &valueTables, size_t index, const std::string &outputFileNamePrefix) {

        const std::string filename = outputFileNamePrefix + (index != 0 ? "#" + to_string(index) : "") + ".eam.fs";

        ofstream eamFsFile(filename);
        if (!eamFsFile.is_open()) {
            JGAP_LOG_ERROR("Could not open " + filename, true);
        }
        eamFsFile << fixed << setprecision(17);

        // Lines 1–3: Comments/metadata.
        eamFsFile << "# UNITS: metal" << std::endl;
        eamFsFile << "# EAM part of a potential tabulated with jGAP"  << std::endl;
        eamFsFile << "# pair_style eam/fs" << std::endl;

        const EamTabulationData& eamTables = valueTables.eamTabulationData.at(index);

        if (eamTables.densityGrids.empty()) {
            JGAP_LOG_AND_THROW("Empty EAM table data");
        }

        // Account for no EAM energy for some element and/or no contribution to density from some element
        // (unlikely to be useful, but being careful is never bad)
        std::set<Species> elementsSet;
        for (const auto& species : eamTables.densityGrids | std::views::keys) {
            elementsSet.insert(species);
        }
        for (const auto& speciesPair : eamTables.eamPairFunctionGrids | std::views::keys) {
            elementsSet.insert(speciesPair.contributor);
            elementsSet.insert(speciesPair.receiver);
        }
        for (const auto& speciesPair : valueTables.pairGrids | std::views::keys) {
            elementsSet.insert(speciesPair.first());
            elementsSet.insert(speciesPair.second());
        }

        // Fix the order, and keep the indexing consistent
        std::vector elements(elementsSet.begin(), elementsSet.end());
        // Line 4: Number of elements (N) followed by each element’s symbol
        eamFsFile << eamTables.densityGrids.size() << " ";
        for (const auto& element: elements) {
            eamFsFile << element << " ";
        }
        eamFsFile << std::endl;

        auto& densityGridSample = eamTables.densityGrids.at(elements[0]);
        for (const auto& grid: eamTables.densityGrids | std::views::values) {
            if (grid.size() != densityGridSample.size()
                || grid.origin != densityGridSample.origin || grid.spacing != densityGridSample.spacing) {
                // Shouldn't happen but just in case
                JGAP_LOG_AND_THROW("Differing EAM density grid setups");
            }
        }

        if (eamTables.eamPairFunctionGrids.empty()) {
            JGAP_LOG_AND_THROW("No EAM pair-function tables");
        }
        auto& pfGridSample = eamTables.eamPairFunctionGrids.begin()->second;
        for (const auto& pfGrid: eamTables.eamPairFunctionGrids | std::views::values) {
            if (pfGrid.size() != pfGridSample.size()
                || pfGrid.origin != pfGridSample.origin || pfGrid.spacing != pfGridSample.spacing) {
                JGAP_LOG_AND_THROW("Differing EAM pair-function grid setups");
            }
        }

        for (const auto& pairEnergyGrid: valueTables.pairGrids | std::views::values) {
            if (pairEnergyGrid.size() != pfGridSample.size()
                || pairEnergyGrid.origin != pfGridSample.origin || pairEnergyGrid.spacing != pfGridSample.spacing) {
                JGAP_LOG_AND_THROW(
                    "EAM pair-function grid setup differs from a pair-energy grid setup"
                    );
            }
        }

        // Line 5: Nrho, drho, Nr, dr, cutoff - so mysterious ..
        eamFsFile << densityGridSample.size() << " ";
        eamFsFile << densityGridSample.spacing << " ";
        eamFsFile << pfGridSample.size() << " ";
        eamFsFile << pfGridSample.spacing << " ";
        eamFsFile << pfGridSample.cutoff() << std::endl;

        // Per-element Sections:
        /*
            * Line: atomic number, mass, lattice constant, lattice type (e.g., fcc, bcc)
            * Embedding function F_\beta(\rho): Nrho tabulated values
            * Density functions \rho_{\alpha\beta}(r): For each element α (total N curves, each with Nr points)
         */
        for (const Species& species1: elements) {
            eamFsFile << static_cast<size_t>(ATOMIC_NUMBERS[species1]) << " ";
            eamFsFile << ATOMIC_MASSES[species1] << " 1.0 ZZZ" << std::endl;

            for (const auto& energyGridSlot : eamTables.getEnergyGridOrNull(species1)) {
                eamFsFile << energyGridSlot.value << std::endl;
            }

            // contributor = species1
            for (const Species& receiver: elements) {
                for (const auto& pairFunctionGridSlot : eamTables.getPairFunctionGridOrNull({species1, receiver})) {
                    eamFsFile << pairFunctionGridSlot.value << std::endl;
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
                    continue;
                }

                for (const auto& pairEnergyGridSlot: valueTables.getOrNull({elements[i], elements[j]})) {
                    eamFsFile << (index == 0 ? pairEnergyGridSlot.value * pairEnergyGridSlot.pos : 0.0) << std::endl;
                }
            }
        }

        eamFsFile.flush();
        eamFsFile.close();
        return filename;
    }

    void TabGapIO::readH5(const std::string &fileName, TabulationData &splineCoefficients) {
        HighFive::File tabGapFile(fileName, HighFive::File::ReadOnly);

        // Read isolated atom energies
        if (tabGapFile.exist("e0")) {
            auto e0Group = tabGapFile.getGroup("e0");
            // Read all attributes except Nelements
            for (const auto &attrName : e0Group.listAttributeNames()) {
                if (attrName == "Nelements") continue;
                double val;
                e0Group.getAttribute(attrName).read(val);
                splineCoefficients.isolatedEnergies[attrName] = val;
            }
        }

        // Read npots to know counts; but we will iterate groups by name
        std::vector<std::string> objectNames = tabGapFile.listObjectNames();
        for (const auto &name : objectNames) {
            if (name == "e0" || name == "npots" || name == "comment1" || name == "comment2") continue;
            // Count number of '-' to distinguish 2b vs 3b
            size_t dashCount = static_cast<size_t>(ranges::count(name, '-'));
            auto group = tabGapFile.getGroup(name);

            if (dashCount == 1) {
                // Pair group: element_i, element_j
                Species iSpec, jSpec;
                group.getAttribute("element_i").read(iSpec);
                group.getAttribute("element_j").read(jSpec);
                SpeciesPair pair{iSpec, jSpec};

                std::vector<double> limits; // origin, cutoff
                tabGapFile.getDataSet(name + std::string("/grid_limits")).read(limits);

                size_t N_originalGrid;
                group.getAttribute("N").read(N_originalGrid);

                std::vector<double> splineCoeffs2b;
                tabGapFile.getDataSet(name + std::string("/energies")).read(splineCoeffs2b);
                if (splineCoeffs2b.size() != N_originalGrid + 2) {
                    JGAP_LOG_AND_THROW("Invalid H5 pair data size for {}-{}", iSpec + std::string("-") + jSpec);
                }

                double originalGridOrigin = limits.at(0);
                double originalGridCutoff = limits.at(1);
                double spacing = (originalGridCutoff - originalGridOrigin) / static_cast<double>(N_originalGrid - 1);
                Grid1d splineGrid(N_originalGrid + 2, spacing, originalGridOrigin - spacing);

                splineGrid.data = std::move(splineCoeffs2b);
                if (splineCoefficients.pairGrids.contains(pair)) {
                    for (const auto& energy: splineCoefficients.pairGrids[pair]) {
                        assert(abs(energy.value) < 1e-20 && "Conflicting 2b energies in files");
                    }
                }
                splineCoefficients.pairGrids[pair] = std::move(splineGrid);

            } else if (dashCount == 2) {
                // Triplet group: element_i, element_j, element_k
                Species iSpec, jSpec, kSpec;
                group.getAttribute("element_i").read(iSpec);
                group.getAttribute("element_j").read(jSpec);
                group.getAttribute("element_k").read(kSpec);

                std::vector<size_t> Nvec; // nR, nR, nAngular
                tabGapFile.getDataSet(name + std::string("/N")).read(Nvec);
                std::vector<double> limits; // rmin,rmin,-1, rmax,rmax,1
                tabGapFile.getDataSet(name + std::string("/grid_limits")).read(limits);
                std::vector<double> splineCoeffs3b;
                tabGapFile.getDataSet(name + std::string("/energies")).read(splineCoeffs3b);

                assert(Nvec[0] == Nvec[1] && "Only symmetric table supported");
                size_t originalGrid_nR = Nvec.at(0);
                size_t originalGrid_nAngular = Nvec.at(2);
                double rminOrigianlGrid = limits.at(0);
                double rmaxOrigianlGrid = limits.at(3);
                double angularMin = limits.at(2); // expected -1
                double angularMax = limits.at(5); // expected 1
                assert(abs(angularMin-(-1)) < 1e-20 && abs(angularMax-1) < 1e-20
                       && "Expected angular limits to be -1 to 1");
                double spacingR = (rmaxOrigianlGrid - rminOrigianlGrid) / static_cast<double>(originalGrid_nR - 1);
                double spacingAngular = (angularMax - angularMin) / static_cast<double>(originalGrid_nAngular - 1);

                assert(splineCoeffs3b.size() == pow(originalGrid_nR + 2, 2) * originalGrid_nAngular
                        && "Number of coefficients mis-matches grid size");

                Grid3d grid(rminOrigianlGrid - spacingR, originalGrid_nR+2, spacingR,
                            angularMin, originalGrid_nAngular+2, spacingAngular);
                grid.dataFlat = std::move(splineCoeffs3b);
                splineCoefficients.tripletGrids[SpeciesTriplet{iSpec, SpeciesPair{jSpec, kSpec}}] = std::move(grid);
            }
        }
    }

    void TabGapIO::readEamFs(const std::string &fileName, TabulationData &splineCoefficients) {
        ifstream inFile(fileName);
        if (!inFile.is_open()) {
            JGAP_LOG_AND_THROW("Could not open {} as .eam.fs", fileName);
        }

        std::string line;
        for (size_t i = 0; i < 3; i++) {
            if (!getLine(inFile, line)) JGAP_LOG_AND_THROW("Invalid EAM/FS: missing comment #{}", i);
        }

        // Elements line: N and symbols
        if (!getLine(inFile, line)) JGAP_LOG_AND_THROW("Invalid EAM/FS: missing elements line");

        std::istringstream iss(line);
        size_t N;
        iss >> N;
        std::vector<Species> elements(N);
        for (size_t i = 0; i < N; i++) iss >> elements[i];
        if (elements.empty()) JGAP_LOG_AND_THROW("Invalid EAM/FS: zero elements");

        // Grid sizes line: Nrho, drho, Nr, dr, cutoff
        if (!getLine(inFile, line)) JGAP_LOG_AND_THROW("Invalid EAM/FS: missing grid spec line");

        size_t Nrho, Nr; double drho, dr, cutoff;
        std::istringstream is2(line);
        is2 >> Nrho >> drho >> Nr >> dr >> cutoff;

        // Prepare a new EAM block in the result
        EamTabulationData eamTables = splineCoefficients.newEamGrids();

        // Per-element sections
        for (const auto& element: elements) {
            // Header line with atomic number, mass, etc. Ignore contents
            if (!getLine(inFile, line)) {
                JGAP_LOG_AND_THROW("Invalid EAM/FS: missing per-element header");
            }

            // Embedding function: Nrho lines
            Grid1d embGrid(Nrho, drho, 0.0);
            for (size_t i = 0; i < Nrho; i++) {
                if (!getLine(inFile, line)) JGAP_LOG_AND_THROW("Invalid EAM/FS: incomplete embedding values");
                embGrid.data[i] = std::stod(line);
            }
            eamTables.densityGrids[element] = BSplineTools::toSplineCoefficients(embGrid);

            // Density functions: for each receiver alpha, Nr points
            for (const auto& receiver: elements) {
                Grid1d rhoGrid(N, dr, 0.0);
                for (size_t i = 0; i < Nr; i++) {
                    if (!getLine(inFile, line)) {
                        JGAP_LOG_AND_THROW("Invalid EAM/FS: incomplete density function");
                    }
                    rhoGrid.data[i] = std::stod(line);
                }
                eamTables.eamPairFunctionGrids[{element, receiver}] = BSplineTools::toSplineCoefficients(rhoGrid);
            }
        }

        // Pair potential tables for i >= j
        for (size_t i = 0; i < N; i++) {
            for (size_t j = 0; j < N; j++) {
                if (i < j) continue;
                Grid1d pairGrid(Nr, dr, 0.0);
                for (size_t k = 0; k < Nr; k++) {
                    if (!getLine(inFile, line)) {
                        JGAP_LOG_AND_THROW("Invalid EAM/FS: incomplete pair potential table");
                    }
                    double phi = std::stod(line);
                    double r = static_cast<double>(k) * dr; // origin assumed 0
                    pairGrid.data[k] = (r > 0.0 ? phi / r : 0.0);
                }
                splineCoefficients.pairGrids[SpeciesPair{elements[i], elements[j]}]
                    = BSplineTools::toSplineCoefficients(pairGrid);
            }
        }
    }
}
#endif
