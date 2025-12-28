#include "core/descriptors/sparsification/HistogramUniformSparsifier.hpp"

#include <utility>

#include "io/log/CurrentLogger.hpp"

namespace jgap {

    HistogramUniformSparsifier::HistogramUniformSparsifier(DataNode params)
            : _kernelParams(std::move(params)) {

        _nSparsePoints = require(_kernelParams, "n_sparse");
        _kernelParams.erase("n_sparse");

        _sparsifiedParamName = require(_kernelParams, "sparse_param");
        _kernelParams.erase("sparse_param");

        if (_kernelParams.contains("seed")) {
            _seed = _kernelParams["seed"];
            _kernelParams.erase("seed");
        } else {
            _seed = 799u;
        }

        if (_kernelParams.contains("grid_dimensions")) {
            _gridDimensions = std::vector<size_t>();
            for (const auto &dim : _kernelParams["grid_dimensions"]) {
                _gridDimensions->push_back(dim.get<size_t>());
            }
            _kernelParams.erase("grid_dimensions");
        }

        if (_kernelParams.contains("min_point")) {
            _minPoint = _kernelParams["min_point"];
            _kernelParams.erase("min_point");
        }
        if (_kernelParams.contains("max_point")) {
            _maxPoint = _kernelParams["max_point"];
            _kernelParams.erase("max_point");
        }
    }

    std::vector<KernelParams> HistogramUniformSparsifier::selectSparsePoints(const std::vector<std::vector<double>> &allPoints) {

        for (const auto &p : allPoints) {
            if (p.size() != allPoints[0].size()) {
                JGAP_LOG_AND_THROW("Sparse points of differing dimensions");
            }
        }

        const std::vector<size_t> gridDimensions = _gridDimensions.value_or(
             std::vector<size_t>(
                    allPoints[0].size(),
                    ceil(pow(_nSparsePoints, 1.0 / static_cast<double>(allPoints[0].size())))
                    )
        );

        if (_minPoint.has_value() && _minPoint->size() != gridDimensions.size()) {
            JGAP_LOG_AND_THROW("Sparsifier's min_point dimensions {} don't match grid dimensions {}",
                                              _minPoint->size(), gridDimensions.size());
        }
        if (_maxPoint.has_value() && _maxPoint->size() != gridDimensions.size()) {
            JGAP_LOG_AND_THROW("Sparsifier's max_point dimensions {} don't match grid dimensions {}",
                                              _maxPoint->size(), gridDimensions.size());
        }

        std::vector<double> minPoint(gridDimensions.size());
        std::vector<double> maxPoint(gridDimensions.size());
        std::vector<double> step(gridDimensions.size());

        for (size_t d = 0; d < gridDimensions.size(); d++) {
            if (_minPoint.has_value()) {
                minPoint[d] = _minPoint.value()[d];
            } else {
                minPoint[d] = std::numeric_limits<double>::max();
                for (const auto &p : allPoints) {
                    minPoint[d] = std::min(minPoint[d], p[d]);
                }
            }
            if (_maxPoint.has_value()) {
                maxPoint[d] = _maxPoint.value()[d];
            } else {
                maxPoint[d] = std::numeric_limits<double>::min();
                for (const auto &p : allPoints) {
                    maxPoint[d] = std::max(maxPoint[d], p[d]+0.0001/*keep all points in bounds*/);
                }
            }
            step[d] = (maxPoint[d] - minPoint[d]) / static_cast<double>(gridDimensions[d]);
        }

        JGAP_LOG_INFO(
                "{}d histogram in range {} - {} with {} long bins:",
                gridDimensions.size(),
                iteratorToString(minPoint.begin(), minPoint.end()),
                iteratorToString(maxPoint.begin(), maxPoint.end()),
                iteratorToString(step.begin(), step.end())
                );

        std::vector<std::vector<double>> sparsePoints;
        std::map<std::vector<size_t>, std::vector<size_t>> usefulGridSlots;

        for (size_t i = 0; i < allPoints.size(); i++) {
            std::vector<size_t> gridSlot{};
            for (size_t d = 0; d < gridDimensions.size(); d++) {
                gridSlot.push_back((allPoints[i][d] - minPoint[d]) / step[d]);
            }

            if (!usefulGridSlots.contains(gridSlot)) {
                usefulGridSlots[gridSlot] = {};
            }
            usefulGridSlots[gridSlot].push_back(i);
        }
        const size_t reps = _nSparsePoints / usefulGridSlots.size();
        const size_t leftover = _nSparsePoints - usefulGridSlots.size() * reps;
        JGAP_LOG_DEBUG(
            "Found {} grid slots containing some points -> attempting to sample each {} times / {} assigned randomly",
            usefulGridSlots.size(), reps, leftover
            );

        std::mt19937 gen(_seed);
        std::vector<std::vector<size_t>> usefulGridSlotsArr = {};
        for (auto &[gridSlot, pointIndices]: usefulGridSlots) {
            usefulGridSlotsArr.push_back(gridSlot);

            std::ranges::shuffle(pointIndices.begin(), pointIndices.end(), gen);
            for (size_t rep = 0; rep < std::min(pointIndices.size(), reps); rep++) {
                sparsePoints.push_back(allPoints[pointIndices[rep]]);
            }
        }

        std::uniform_int_distribution<> indexDist(0, usefulGridSlotsArr.size() - 1);

        while (sparsePoints.size() < _nSparsePoints) {
            const std::vector<size_t>& gridSlot = usefulGridSlotsArr[indexDist(gen)];

            std::vector<double> point(gridDimensions.size());
            for (size_t d = 0; d < gridDimensions.size(); d++) {
                std::uniform_real_distribution<> marginDist(0, step[d]);
                point[d] = minPoint[d] + step[d] * static_cast<double>(gridSlot[d]) + marginDist(gen);
            }

            sparsePoints.push_back(point);
            JGAP_LOG_DEBUG(iteratorToString(point.begin(), point.end()));
        }

        std::vector<KernelParams> sparseKernelsParams;
        for (const auto& point: sparsePoints) {
            sparseKernelsParams.push_back(_kernelParams);
            if (point.size() == 1) {
                sparseKernelsParams.back()[_sparsifiedParamName] = point[0];
            } else {
                sparseKernelsParams.back()[_sparsifiedParamName] = point;
            }
        }

        return sparseKernelsParams;
    }
}
