#include "app/App.hpp"

#include <Eigen/Dense>
#include <Version.hpp>
#include <core/fit/Fit.hpp>
#include <core/potentials/Potential.hpp>
#include <cstdlib>
#include <filesystem>
#include <io/log/CurrentLogger.hpp>
#include <queue>
#include <string>
#include <thread>
#include <vector>

#include "app/FitApp.hpp"

using namespace std;

namespace jgap {
    int App::main(int argc, char** argv) {

        CurrentLogger::get()->info("jGAP app v{}", JGAP_VERSION);

        Eigen::setNbThreads(thread::hardware_concurrency());

        vector<string> args{};
        vector<string> flags{};
        for (size_t i = 1; i < argc; i++) {
            if (argv[i][0] == '-') {
                flags.emplace_back(argv[i]);
            } else {
                args.emplace_back(argv[i]);
            }
        }

        if (args.empty() && flags.empty()) {
            CurrentLogger::get()->error("No input | use --help to see available options");
            return EXIT_FAILURE;
        }

        optional<vector<AtomicStructure>> contextXyz;
        optional<TabulationParams> contextTabulationParams;
        shared_ptr<Fit> contextFit;
        shared_ptr<Potential> contextPotential;

        while (!args.empty() || !flags.empty()) {

        }

        if (!flags.empty()) {
            if (flags[0] == "-h" || flags[0] == "-help" || flags[0] == "--help") {
                printHelp();
                return EXIT_SUCCESS;
            }

            if (flags[0] == "-v" || flags[0] == "-version" || flags[0] == "--version") {
                CurrentLogger::get()->info("jGAP app v{}", JGAP_VERSION);
                return EXIT_SUCCESS;
            }
        } else {
            // Auto-detect based on file extension

            if (args.size() == 1) {

                if (args[0].ends_with(".json")) {
                    return EXIT_SUCCESS;
                }


            } else {
            }
        }



        return EXIT_SUCCESS;
    }

    void App::printHelp() {

    }
    shared_ptr<Potential> App::autoDetectPotential() {}
    optional<nlohmann::json> App::autoDetectParams() {}
    optional<vector<AtomicStructure>> App::autoDetectXyz() {}
} // namespace jgap
