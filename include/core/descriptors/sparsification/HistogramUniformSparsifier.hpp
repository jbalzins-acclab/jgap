#ifndef HISTOGRAMUNIFORMSPARSIFIER_HPP
#define HISTOGRAMUNIFORMSPARSIFIER_HPP

#include <nlohmann/json.hpp>
#include <cmath>
#include <random>
#include <set>

#include "Sparsifier.hpp"
#include "io/parse/ParserRegistry.hpp"
#include "utils/Utils.hpp"

namespace jgap {

    class HistogramUniformSparsifier : public Sparsifier {
    public:
        HistogramUniformSparsifier(nlohmann::json params);
        vector<nlohmann::json> selectSparsePoints(const vector<vector<double>> &allPoints) override;
    private:
        size_t _seed;
        size_t _nSparsePoints;
        string _sparsifiedParamName;
        nlohmann::json _kernelParams;
        optional<vector<size_t>> _gridDimensions;
        optional<vector<double>> _minPoint;
        optional<vector<double>> _maxPoint;
    };

    REGISTER_PARSER("histogram_uniform", Sparsifier, HistogramUniformSparsifier);
}

#endif
