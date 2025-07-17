#ifndef TABULATE_HPP
#define TABULATE_HPP

#include <memory>
#include <nlohmann/json.hpp>

#include "core/potentials/Potential.hpp"

using namespace std;

namespace jgap {
    class Tabulate {
    public:
        static void tabulate(const shared_ptr<Potential> &potential,
                             const nlohmann::json& params,
                             const string &outputFileNamePrefix);
    private:
        static TabulationParams parse(const nlohmann::json& params);
        static vector<double> make2bGrid(const nlohmann::json& params);
        static vector<vector<vector<Vector3>>> make3bGrid(const nlohmann::json& params);
        static vector<double> toSplineCoefficients(const vector<double>& energies, double spacing);
        static vector<double> toSplineCoefficients(const vector<vector<vector<double>>>& energies, const Vector3 &spacing);
        static void writeH5(const shared_ptr<Potential> &potential,
                            const nlohmann::json& params,
                            const TabulationParams &tabulationParams,
                            const TabulationData &tabulationData,
                            const string &outputFileNamePrefix);
        static void writeEamFs(const shared_ptr<Potential> &potential,
                               const nlohmann::json &params,
                               const TabulationParams &tabulationParams,
                               const EamTabulationData& data,
                               const string &outputFileNamePrefix,
                               const map<SpeciesPair, vector<double>> &pairEnergies,
                               bool writePairEnergies);
    };
}

#endif //TABULATE_HPP
