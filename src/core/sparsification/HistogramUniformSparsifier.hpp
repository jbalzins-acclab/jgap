#ifndef HISTOGRAMUNIFORMSPARSIFIER_HPP
#define HISTOGRAMUNIFORMSPARSIFIER_HPP

#include <cmath>
#include <random>
#include <set>

#include "Sparsifier.hpp"
#include "io/parse/ParserRegistry.hpp"
#include "utils/Utils.hpp"

namespace jgap {

    class HistogramUniformSparsifier : public Sparsifier {
    public:
        SETUP_PARSER(Sparsifier, HistogramUniformSparsifier, histogram_uniform)

        HistogramUniformSparsifier();

        std::vector<KernelParams> selectSparsePoints(const std::vector<std::vector<double>> &allPoints) override;

    private:
        size_t _seed;
        size_t _nSparsePoints;
        std::string _sparsifiedParamName;
        DataNode _kernelParams;
        std::optional<std::vector<size_t>> _gridDimensions;
        std::optional<std::vector<double>> _minPoint;
        std::optional<std::vector<double>> _maxPoint;
    };
}

#endif
