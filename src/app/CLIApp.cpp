#include "app/CLIApp.hpp"

#include <Eigen/Dense>
#include <Version.hpp>
#include <core/fit/Fit.hpp>
#include <core/potentials/Potential.hpp>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <io/log/CurrentLogger.hpp>
#include <queue>
#include <string>
#include <thread>
#include <vector>
#include <oneapi/tbb/parallel_for.h>

#include "io/convert/QuipXmlConverter.hpp"

namespace jgap {

    int CLIApp::main(int argc, char** argv) {

#ifdef SET_LOGGER
        CurrentLogger::setLogger(SET_LOGGER);
#endif

        JGAP_LOG_INFO("jGAP app v{}", JGAP_VERSION);

        Eigen::setNbThreads(std::thread::hardware_concurrency());

        if (argc < 2) {
            JGAP_LOG_ERROR("No arguments provided");
            return EXIT_FAILURE;
        }

        if (HELP_FLAGS.contains(argv[1])) {
            printHelp();
            return EXIT_SUCCESS;
        }

        if (VERSION_FLAGS.contains(argv[1])) {
            // always in first message
            return EXIT_SUCCESS;
        }

        if (processedAsConversion(argc, argv)) return EXIT_SUCCESS;
        if (processedAsFit(argc, argv)) return EXIT_SUCCESS;
        if (processedAsTabulation(argc, argv)) return EXIT_SUCCESS;
        if (processedAsActionJson(argc, argv)) return EXIT_SUCCESS;
        if (processedAsPrediction(argc, argv)) return EXIT_SUCCESS;
        if (processedAsRmseTesting(argc, argv)) return EXIT_SUCCESS;

        JGAP_LOG_ERROR("Failed to process provided arguments | use --help to compare with available options");
        return EXIT_FAILURE;
    }

    bool CLIApp::processedAsConversion(int argc, char **argv) {

        std::string quipPotFileName;
        if (CONVERT_FLAGS.contains(argv[1])) {
            if (argc < 3) {
                JGAP_LOG_AND_THROW("File to be converted not specified");
            }
            quipPotFileName = argv[2];
        } else if (std::string(argv[1]).ends_with(".xml")) {
            quipPotFileName = argv[1];
        } else {
            return false;
        }

        auto potential = convertQuipXml(quipPotFileName);

        std::optional<std::string> outputFileName{};
        if (argv[argc - 1] != quipPotFileName) {
            outputFileName = argv[argc - 1];
        } else if (quipPotFileName.ends_with(".xml")) {
            outputFileName = quipPotFileName.substr(0, quipPotFileName.size() - 4) + ".jgap.json";
        }

        savePotential(potential, outputFileName);

        return true;
    }

    bool CLIApp::processedAsFit(int argc, char **argv) {
        if (FIT_FLAGS.contains(argv[1])) return false;

        if (argc < 3) {
            JGAP_LOG_AND_THROW("Not enough arguments for a fit");
        }

        std::ifstream paramFile(argv[2]);
        if (!paramFile.is_open()) {
            JGAP_LOG_AND_THROW("Cannot open param file {}", argv[2]);
        }
        nlohmann::json fitParams;
        try {
            paramFile >> fitParams;
        } catch (const std::exception& e) {
            JGAP_LOG_AND_THROW("{} not a valid json: {}", argv[2], e.what());
        }

        JGAP_LOG_INFO("Fitting potential");
        auto potential = fit(fitParams);

        std::optional<std::string> fitOutputFileName{};
        optionallySet(fitOutputFileName, fitParams, "out");

        JGAP_LOG_INFO("Saving potential");
        savePotential(potential, fitOutputFileName);

        if (fitParams.contains("tabulation")) {
            nlohmann::json tabulationParams = fitParams["tabulation"];

            tabulationParams["type"] = tabulationParams.value("type", SimpleTabulation::TYPE);

            auto tabulation = REGISTRY_GET(Tabulation, fitParams["tabulation"]);

            JGAP_LOG_INFO("Tabulating potential");
            auto tabulatedPotential = tabulation->tabulate(potential);

            std::optional<std::string> tabOutputFileName{};
            optionallySet(tabOutputFileName, tabulationParams, "out");
            if (!tabOutputFileName.has_value() && fitOutputFileName.has_value()) {
                tabOutputFileName = "tab-" + *fitOutputFileName;
            }

            JGAP_LOG_INFO("Saving markup for the tabulated version");
            savePotential(potential, tabOutputFileName);
        }

        return true;
    }

    bool CLIApp::processedAsRmseTesting(int argc, char **argv) {
        if ()
    }

    std::shared_ptr<Potential> CLIApp::fit(nlohmann::json& fitParamNode) {

        // ------------------------ READ PARAMS AND PREPARE -------------------------------
        fitParamNode["type"] = fitParamNode.value("type", "composite");
        JGAP_LOG_DEBUG("{} fit", fitParamNode["type"].get<std::string>());
        auto fit = REGISTRY_GET(Fit, fitParamNode);

        JGAP_LOG_DEBUG(
            "Reading training data from {}", require(fitParamNode, "training_data_xyz").get<std::string>()
            );
        auto trainingData = readXyz(require(fitParamNode, "training_data_xyz"));

        // ------------------------ FIT -------------------------------

        JGAP_LOG_INFO("Fitting {} potential", fitParamNode["type"].get<std::string>());
        JGAP_LOG_DEBUG("Fit params: {}", fitParamNode.dump());

        std::shared_ptr<Potential> resultingPotential = fit->fit(trainingData);
        JGAP_LOG_INFO("Fit done");

        return resultingPotential;


        // ------------------------ TEST -------------------------------

        if (fitParams.value("do_test", true)) {
            JGAP_LOG_INFO("Running RMSE tests");

            std::set<std::string> groupTestDataBy = fitParams.value(
                "group_test_data_by",
                std::set<std::string>{"overall", "config_type"}
                );

            std::set<std::string> testFiles = fitParams.value("test_files", std::set<std::string>{});
            testFiles.insert(fitParams["training_data_xyz"]);

            for (const std::string& testFile: testFiles) {
                JGAP_LOG_INFO("Testing {}", testFile);
                auto testingData = readXyz(testFile, resultingPotential->getCutoff().maxOverall());

                std::vector<Predictions> predictions(testingData.size());
                tbb::parallel_for(0uz, testingData.size(), [&](const size_t i) {
                    predictions[i] = resultingPotential->predict(testingData[i]);
                });

                for (const std::string& propertyName: groupTestDataBy) {
                    std::map<std::string, std::vector<double> > energyDifferences;
                    std::map<std::string, std::vector<double> > forceDifferences;
                    std::map<std::string, std::vector<std::array<Vector3, 3> > > virialDifferences;

                    for (size_t i = 0; i < testingData.size(); i++) {

                        std::string propertyValue = GET_OR_DEFAULT(testingData[i].properties, propertyName, "-");
                        if (!energyDifferences.contains(propertyValue)) {
                            energyDifferences[propertyValue] = {};
                            forceDifferences[propertyValue] = {};
                            virialDifferences[propertyValue] = {};
                        }

                        double nAtoms = testingData[i].size();
                        if (predictions[i].energy.has_value() && testingData[i].energy.has_value()) {
                            energyDifferences[propertyValue].push_back(
                                (predictions[i].energy.value() - testingData[i].energy.value()) / nAtoms
                                );
                        }

                        if (predictions[i].virials.has_value() && testingData[i].virials.has_value()) {
                            virialDifferences[propertyValue].push_back({
                                (predictions[i].virials.value()[0] - testingData[i].virials.value()[0]) / nAtoms,
                                (predictions[i].virials.value()[1] - testingData[i].virials.value()[1]) / nAtoms,
                                (predictions[i].virials.value()[2] - testingData[i].virials.value()[2]) / nAtoms
                            });
                        }

                        if (!testingData[i].forces.has_value() || !predictions[i].forces.has_value()) continue;
                        for (size_t j = 0; j < testingData[i].size(); j++) {
                            forceDifferences[propertyValue].push_back(
                                ((*predictions[i].forces)[j] - (*testingData[i].forces)[j]).len()
                            );
                        }
                    }

                    printErrors(
                        outputFileName, testFile, propertyName,
                        energyDifferences, forceDifferences, virialDifferences
                        );
                }
            }
        }
    }

    void CLIApp::savePotential(const std::shared_ptr<Potential> &potential, std::optional<std::string> outputFileName) {

        if (!outputFileName.has_value()) {
            outputFileName = potential->getType() + potential->hashString() + ".jgap.json";
            JGAP_LOG_WARN(
                "Output file name not specified for a potential. Using default value {} to save",
                outputFileName.value()
                );
        }

        auto output = potential->serialize();
        output["type"] = potential->getType();

        JGAP_LOG_INFO("Saving potential to: {}", *outputFileName);
        std::ofstream outputFile(*outputFileName);
        if (!outputFile.is_open()) {
            JGAP_LOG_AND_THROW("Could not open output file {}", *outputFileName);
        }
        outputFile << output;
        outputFile.flush();
        outputFile.close();
    }

    std::shared_ptr<Potential> CLIApp::convertQuipXml(std::string quipFileName) {
        JGAP_LOG_INFO("Converting QUIP potential {} to jGAP", quipFileName);

        pugi::xml_document quipDocument;
        if (!quipDocument.load_file(quipFileName.c_str())) {
            JGAP_LOG_AND_THROW("Cannot open the input quip file: {}", quipFileName);
        }

        auto result = QuipXmlConverter::transform(quipDocument.document_element());

        JGAP_LOG_INFO("Conversion complete");
        return result;
    }

    void CLIApp::printHelp() {

    }

    void printErrors(const std::string& potFileName, const std::string &testFile, const std::string &groupedBy,
                     std::map<std::string, std::vector<double> > &energyDifferences,
                     std::map<std::string, std::vector<double> > &forceDifferences,
                     std::map<std::string, std::vector<std::array<Vector3, 3> > > &virialDifferences) {

        JGAP_LOG_DEBUG("Grouping {} test {}'s predictions by {}", testFile, potFileName, groupedBy);
        std::ofstream reportFile(format("{}-{}-errors.{}.csv", potFileName, split(testFile, '/').back(), groupedBy));

        reportFile << "groupedBy,"
                      "N_energies,E_RMSE(meV/atom),E_STDE(meV),"
                      "N_Forces,F_RMSE(eV/A),F_STDE(eV/A),"
                      "N_virials,"
                      "Virial_isotropic_RMSE(eV/atom),"
                      "Virial_shear_RMSE(eV/atom)"
                   << std::endl;

        for (const std::string& propertyValue: energyDifferences | std::views::keys) {
            reportFile << propertyValue << ",";

            reportFile << energyDifferences[propertyValue].size() << ",";
            reportFile << rms(energyDifferences[propertyValue]) * 1e3 << ",";
            reportFile << deviation(energyDifferences[propertyValue]) * 1e3 << ",";

            reportFile << forceDifferences[propertyValue].size() << ",";
            reportFile << rms(forceDifferences[propertyValue]) << ",";
            reportFile << deviation(forceDifferences[propertyValue]) << ",";

            std::vector<double> isotropicVirialErr, anisotropicVirialErr;
            for (const auto& virialErr: virialDifferences[propertyValue]) {
                isotropicVirialErr.push_back((virialErr[0].x + virialErr[1].x + virialErr[2].x) / 3.0);
                anisotropicVirialErr.push_back(
                    (virialErr[0].y + virialErr[0].z + virialErr[1].x +
                        virialErr[1].z + virialErr[2].x + virialErr[2].y) / 6.0
                    );
            }
            reportFile << virialDifferences[propertyValue].size() << ",";
            reportFile << rms(isotropicVirialErr) << ",";
            reportFile << rms(anisotropicVirialErr);

            reportFile << std::endl;
        }

        reportFile.flush();
        reportFile.close();
    }
}
