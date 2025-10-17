#include <exception>
#include <iostream>
#include <ostream>

#include "app/App.hpp"

int main(int argc, char** argv) {
#ifndef DEBUG
    try {
        return jgap::App::main(argc, argv);
    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }
#else
    // Not to lose the exception when doing a debug run.
    return jgap::App::main(argc, argv);
#endif
}