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
#include <yaml-cpp/yaml.h>

#include "app/ProcessingPipeline.hpp"
#include "app/instructions/LoadInstruction.hpp"
#include "io/convert/QuipXmlConverter.hpp"
#include "io/convert/YamlConverter.hpp"

namespace jgap {
    int CLIApp::main(int argc, char **argv) {
#ifdef SET_LOGGER
        CurrentLogger::setLogger(SET_LOGGER);
#endif

        JGAP_LOG_INFO("jGAP app v{}", JGAP_VERSION);

        Eigen::setNbThreads(std::thread::hardware_concurrency());

        if (argc < 2) {
            JGAP_LOG_ERROR("No arguments provided");
            return EXIT_FAILURE;
        }

        if (HELP_FLAGS_.contains(argv[1])) {
            help();
            return EXIT_SUCCESS;
        }

        if (VERSION_FLAGS_.contains(argv[1])) {
            // always in the first log info
            return EXIT_SUCCESS;
        }

        if (PRINT_INSTRUCTIONS_.contains(argv[1])) {

        }

        if (std::string(argv[1]).ends_with(".xml")) {
            return processAsConversion(argc, argv);
        }

        if (std::string(argv[1]).ends_with(".yaml")) {
            return processAsConversion(argc, argv);
        }

        return processAsInlineInstructions(argc, argv);
    }

    void CLIApp::help() {

    }

    void CLIApp::printInstructions() {

    }

    int CLIApp::processAsConversion(int argc, char **argv) {
        JGAP_LOG_INFO("First argument is an xml file => processing as conversion");
        if (argc > 3) {
            JGAP_LOG_ERROR("Too many arguments for conversion | use to_be_converted.xml output.yaml");
            return EXIT_FAILURE;
        }
        std::string quip_pot_file_name = argv[1];

        std::optional<std::string> output_file_name{};
        if (argv[argc - 1] != quip_pot_file_name) {
            output_file_name = argv[argc - 1];
        } else if (quip_pot_file_name.ends_with(".xml")) {
            output_file_name = quip_pot_file_name.substr(0, quip_pot_file_name.size() - 4) + ".jgap.json";
        }

        ProcessingPipeline pipeline(std::vector<std::shared_ptr<PipelineInstruction>>{
            std::make_shared<LoadInstruction>(quip_pot_file_name, ""),

        });
        auto potential = convertQuipXml(quip_pot_file_name);
        savePotential(potential, output_file_name);

        return true;
    }

    int CLIApp::processAsInlineInstructions(int argc, char **argv) {
        JGAP_LOG_INFO("Processing inline instructions");
        std::string instructionsString;
        for (int i = 1; i < argc; i++) {
            instructionsString += argv[i];
            instructionsString += " ";
        }
        return processInstructions(instructionsString);
    }

    int CLIApp::processAsInFileInstructions(int argc, char **argv) {
        JGAP_LOG_INFO("Processing {} as instructions file", argv[1]);
        if (argc != 2) {
            JGAP_LOG_ERROR("Expected exactly one argument - instructions file");
            return EXIT_FAILURE;
        }
        std::ifstream inFile(argv[1]);
        if (!inFile.is_open()) {
            JGAP_LOG_ERROR("Could not open file {}", argv[1]);
        }

        std::ostringstream instructionsStringStream;
        instructionsStringStream << inFile.rdbuf();

        return processInstructions(instructionsStringStream.str());
    }

    int CLIApp::processInstructions(const std::string &instructionsString) {

        YAML::Node instructionsYamlNode;
        try {
            instructionsYamlNode = YAML::Load(instructionsString);
        } catch (YAML::ParserException& e) {
            JGAP_LOG_ERROR("Could not parse yaml instructions: {}", e.what());
            return EXIT_FAILURE;
        }

        DataNode instructionsNode = YamlConverter::fromYaml(instructionsYamlNode);

        if (instructionsNode.type == DataNode::Type::ARRAY) {

        } else if (instructionsNode.type == DataNode::Type::OBJECT) {

        } else if (instructionsNode.type == DataNode::Type::STRING) {

        } else {
            JGAP_LOG_ERROR();
            return EXIT_FAILURE;
        }

        return EXIT_SUCCESS;
    }
}
