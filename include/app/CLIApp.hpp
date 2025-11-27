#ifndef JGAP_CLIAPP_HPP
#define JGAP_CLIAPP_HPP

#include <core/potentials/Potential.hpp>
#include <optional>

namespace jgap {
    class CLIApp {
    public:
        static int main(int argc, char** argv);

    private:
        inline static const std::set<std::string> VERSION_FLAGS{"-v", "-version", "--version"};
        inline static const std::set<std::string> HELP_FLAGS{"-h", "-help", "--help"};
        inline static const std::set<std::string> CONVERT_FLAGS{"-c", "-convert", "--convert"};
        inline static const std::set<std::string> FIT_FLAGS{"-f", "-fit", "--fit"};
        inline static const std::set<std::string> TABULATE_FLAGS{"-t", "-tab", "--tab", "-tabulate", "--tabulate"};
        inline static const std::set<std::string> RMS_TEST_FLAGS{"-rmse_test", "--rmse_test"};

        // jgap (-convert) quip_pot.xml (conv_pot.jgap.json)
        static bool processedAsConversion(int argc, char** argv);
        // jgap -fit params.json (out_pot_name.jgap.json) (<- tabulate if contains "tabulation"
        // -> if either node has no output pattern spec => deduce from where it is specified /
        // if neither -> save to )
        // TODO: jgap -fit training_data.xyz (out_pot_name.jgap.json) => auto-test on training data
        static bool processedAsFit(int argc, char** argv);
        // jgap -tab params.json
        // jgap -tab pot.jgap.json - sensitive to .jgap before .json
        static bool processedAsTabulation(int argc, char** argv);
        // jgap -predict pot.jgap.json to_be_predicted.xyz write_res_to.xyz
        static bool processedAsPrediction(int argc, char** argv);
        // jgap -rmse_test pot.jgap.json  (-group_by=config_type,header_name2,...) test.xyz files.xyz ...
        static bool processedAsRmseTesting(int argc, char** argv);
        // jgap do_stuff.json
        /*
         * {
         * ("fit":{...} or "fits":[{train_data=,...},{(last_fit_as_external=false/true),...},{...}]),
         * ("tabulation":{(potential_file or pot) or lastest fit result if spec, ...})
         * }
         */

        static std::shared_ptr<Potential> fit(nlohmann::json& fitParamNode);
        static void doTest(std::shared_ptr<Potential> potential, );
        static void printErrors(const std::string& potFileName, const std::string& testFile, const std::string& groupedBy,
                                std::map<std::string, std::vector<double>>& energyDifferences,
                                std::map<std::string, std::vector<double>>& forceDifferences,
                                std::map<std::string, std::vector<std::array<Vector3, 3>>>& virialDifferences);

        static std::shared_ptr<Potential> convertQuipXml(std::string quipFileName);
        static void savePotential(const std::shared_ptr<Potential> &potential, std::optional<std::string> outputFileName);

        static void printHelp();
    };
}

#endif