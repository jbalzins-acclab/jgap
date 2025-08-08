#ifndef SAMPLESPACEUNIFORMSPARSIFIER_HPP
#define SAMPLESPACEUNIFORMSPARSIFIER_HPP

#include <nlohmann/json.hpp>
#include <cmath>
#include <random>
#include <set>

#include "Sparsifier.hpp"
#include "io/parse/ParserRegistry.hpp"
#include "utils/Utils.hpp"

namespace jgap {

    class SampleSpaceUniformSparsifier : public Sparsifier {
    public:
        explicit SampleSpaceUniformSparsifier(const nlohmann::json& params);
        vector<vector<double>> selectSparsePoints(const vector<vector<double>> &allPoints) override;
    private:
        size_t _nSparsePoints;
        size_t _seed;
        optional<vector<size_t>> _gridDimensions;
    };

    REGISTER_PARSER("sample_space_uniform", Sparsifier, SampleSpaceUniformSparsifier);

}

#endif
