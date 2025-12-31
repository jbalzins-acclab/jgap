#ifndef JGAP_RMSEINSTRUCTION_HPP
#define JGAP_RMSEINSTRUCTION_HPP

#include "PipelineInstruction.hpp"

namespace jgap {
    class RmseInstruction : public PipelineInstruction {
    public:
        RmseInstruction();
        void execute(PipelineState &state) override;
    };
}
    void printErrors(const std::string &potFileName, const std::string &testFile, const std::string &groupedBy,
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

        for (const std::string &propertyValue: energyDifferences | std::views::keys) {
            reportFile << propertyValue << ",";

            reportFile << energyDifferences[propertyValue].size() << ",";
            reportFile << rms(energyDifferences[propertyValue]) * 1e3 << ",";
            reportFile << deviation(energyDifferences[propertyValue]) * 1e3 << ",";

            reportFile << forceDifferences[propertyValue].size() << ",";
            reportFile << rms(forceDifferences[propertyValue]) << ",";
            reportFile << deviation(forceDifferences[propertyValue]) << ",";

            std::vector<double> isotropicVirialErr, anisotropicVirialErr;
            for (const auto &virialErr: virialDifferences[propertyValue]) {
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

if (fitParams.value("do_test", true)) {
            JGAP_LOG_INFO("Running RMSE tests");

            std::set<std::string> groupTestDataBy = fitParams.value(
                "group_test_data_by",
                std::set<std::string>{"overall", "config_type"}
            );

            std::set<std::string> testFiles = fitParams.value("test_files", std::set<std::string>{});
            testFiles.insert(fitParams["training_data_xyz"]);

            for (const std::string &testFile: testFiles) {
                JGAP_LOG_INFO("Testing {}", testFile);
                auto testingData = readXyz(testFile, resultingPotential->getCutoff().maxOverall());

                std::vector<Predictions> predictions(testingData.size());
                tbb::parallel_for(0uz, testingData.size(), [&](const size_t i) {
                    predictions[i] = resultingPotential->predict(testingData[i]);
                });

                for (const std::string &propertyName: groupTestDataBy) {
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

#endif