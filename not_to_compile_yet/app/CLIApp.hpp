#ifndef JGAP_CLIAPP_HPP
#define JGAP_CLIAPP_HPP

#include <core/potentials/Potential.hpp>
#include <optional>

namespace jgap {
    class CLIApp {
    public:
        static int main(int argc, char **argv);

    private:
        inline static const std::set<std::string> VERSION_FLAGS_{"-v", "-version", "--version"};
        inline static const std::set<std::string> HELP_FLAGS_{"-h", "-help", "--help"};
        inline static const std::set<std::string> PRINT_INSTRUCTIONS_{
            "-pi", "-print_instructions", "--print_instructions"
        };

        static void help();

        static void printInstructions();

        static int processAsConversion(int argc, char **argv);

        // predict = jgap [load pot.yaml, load: data.xyz, predict]
        static int processAsInlineInstructions(int argc, char **argv);
        static int processAsInFileInstructions(int argc, char **argv);

        static int processInstructions(const std::string& instructionsString);
    };
}

#endif