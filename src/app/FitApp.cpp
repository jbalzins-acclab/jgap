#include <core/fit/QRGapFit.hpp>
#include <utils/Utils.hpp>
#include <fstream>
#include <nlohmann/json.hpp>
#include <thread>
#include <Eigen/Dense>
#include <ParserRegistryAuto.hpp>
#include <tbb/parallel_for_each.h>

#include <Version.hpp>

using namespace std;

namespace jgap {

        void fit() {

        // ------------------------ READ PARAMS AND PREPARE -------------------------------

            string paramFileName = argv[1];
            ifstream paramFile(paramFileName);
            if (!paramFile.is_open()) {
                jgap::CurrentLogger::get()->error("Cannot open param file {}", paramFileName);
                return EXIT_FAILURE;
            }
            nlohmann::json fitParams;
            paramFile >> fitParams;

            string outputFileName = fitParams["output_file"];
            jgap::CurrentLogger::get()->info("Output file name: {}", outputFileName);

            jgap::CurrentLogger::get()->info("Picking fit logic");
            fitParams["type"] = fitParams.value("type", "composite");
            auto fit = jgap::ParserRegistry<jgap::Fit>::get(fitParams);

            jgap::CurrentLogger::get()->info("Reading training data");
            auto trainingData = jgap::readXyz(fitParams["training_data_xyz"]);

            if (!fitParams.value("use_virials", true)) {
                jgap::CurrentLogger::get()->warn("Not using virials in training data");

                for (auto& structure: trainingData) {
                    structure.virials.reset();
                }
            }

            // ------------------------ FIT -------------------------------

            jgap::CurrentLogger::get()->info(
                "Fitting {} potential with params from file {}: {}",
                fitParams["type"].dump(), paramFileName, fitParams.dump()
                );
            shared_ptr<jgap::Potential> resultingPotential = fit->fit(trainingData);
            jgap::CurrentLogger::get()->info("Main fit done");

            // ------------------------ SAVE -------------------------------

            jgap::CurrentLogger::get()->info("Saving resulting potential data to {}", outputFileName);
            auto output = resultingPotential->serialize();
            output["type"] = resultingPotential->getType();

            ofstream outFile(outputFileName);
            outFile << output.dump(4) << endl;
            outFile.flush();
            outFile.close();

            jgap::CurrentLogger::get()->info("Fit complete");

            if (fitParams.value("do_test", true)) {
                jgap::CurrentLogger::get()->info("Running RMSE tests");

                set<string> groupTestDataBy = fitParams.value(
                    "group_test_data_by",
                    set<string>{"overall", "config_type"}
                    );

                set<string> testFiles = fitParams.value("test_files", set<string>{});
                testFiles.insert(fitParams["training_data_xyz"]);

                for (const string& testFile: testFiles) {
                    jgap::CurrentLogger::get()->info("Testing {}", testFile);
                    auto testingData = jgap::readXyz(testFile, resultingPotential->getCutoff());

                    vector<jgap::PotentialPrediction> predictions(testingData.size());
                    tbb::parallel_for(0uz, testingData.size(), [&](const size_t i) {
                        predictions[i] = resultingPotential->predict(testingData[i]);
                    });

                    for (const string& propertyName: groupTestDataBy) {
                        map<string, vector<double> > energyDifferences;
                        map<string, vector<double> > forceDifferences;
                        map<string, vector<array<jgap::Vector3, 3> > > virialDifferences;

                        for (size_t i = 0; i < testingData.size(); i++) {

                            string propertyValue = GET_OR_DEFAULT(testingData[i].properties, propertyName, "-");
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

void printErrors(const string& potFileName, const string &testFile, const string &groupedBy,
                 map<string, vector<double> > &energyDifferences,
                 map<string, vector<double> > &forceDifferences,
                 map<string, vector<array<jgap::Vector3, 3> > > &virialDifferences) {

    jgap::CurrentLogger::get()->info("Testing {}", format("{}-{}-errors.{}.csv",
        potFileName,
        jgap::split(testFile, '/').back(),
        groupedBy
        ));
    ofstream reportFile(format("{}-{}-errors.{}.csv",
        potFileName,
        jgap::split(testFile, '/').back(),
        groupedBy
        ));
    reportFile << "groupedBy,"
                  "N_energies,E_RMSE(meV/atom),E_STDE(meV),"
                  "N_Forces,F_RMSE(eV/A),F_STDE(eV/A),"
                  "N_virials,"
                  "Virial_isotropic_RMSE(eV/atom),"
                  "Virial_shear_RMSE(eV/atom)" << endl;

    for (const string& propertyValue: energyDifferences | views::keys) {
        reportFile << propertyValue << ",";

        reportFile << energyDifferences[propertyValue].size() << ",";
        reportFile << rms(energyDifferences[propertyValue])*1e3 << ",";
        reportFile << deviation(energyDifferences[propertyValue])*1e3 << ",";

        reportFile << forceDifferences[propertyValue].size() << ",";
        reportFile << rms(forceDifferences[propertyValue]) << ",";
        reportFile << deviation(forceDifferences[propertyValue]) << ",";

        vector<double> isotropicVirialErr, anisotropicVirialErr;
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

        reportFile << endl;
    }

    reportFile.flush();
    reportFile.close();
}

double rms(const vector<double> &x) {
    if (x.empty()) return 0.0;
    double res = 0.0;
    for (double i : x) {
        res += i*i;
    }
    res = sqrt(res / static_cast<double>(x.size()));
    return res;
}

double deviation(const vector<double> &x) {
    if (x.empty()) return 0.0;
    double mean = 0.0;
    for (double i : x) {
        mean += i;
    }
    mean /= static_cast<double>(x.size());
    double variance = 0.0;
    for (double i : x) {
        variance += (i - mean) * (i - mean);
    }
    variance = variance / static_cast<double>(x.size());
    return sqrt(variance);
}
