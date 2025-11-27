#include <exception>
#include <iostream>
#include <ostream>

#include "app/CLIApp.hpp"

int main(int argc, char** argv) {
#ifndef DEBUG
    try {
        return jgap::CLIApp::main(argc, argv);
    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }
#else
    // Not to lose the exception when doing a debug run.
    return jgap::CLIApp::main(argc, argv);
#endif
}