#ifndef JGAP_APP_HPP
#define JGAP_APP_HPP

#include <core/potentials/Potential.hpp>
#include <optional>

using namespace std;

namespace jgap {
    class App {
    public:
        static int main(int argc, char** argv);
    private:
        static void printHelp();
        static shared_ptr<Potential> autoDetectPotential();
        static optional<nlohmann::json> autoDetectParams();
        static optional<vector<AtomicStructure>> autoDetectXyz();
    };
}

#endif
