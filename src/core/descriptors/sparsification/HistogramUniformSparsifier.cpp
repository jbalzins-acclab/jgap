#include "core/descriptors/sparsification/HistogramUniformSparsifier.hpp"

#include <utility>

namespace jgap {

    HistogramUniformSparsifier::HistogramUniformSparsifier(nlohmann::json params)
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
            _gridDimensions = vector<size_t>();
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

    vector<KernelParams> HistogramUniformSparsifier::selectSparsePoints(const vector<vector<double>> &allPoints) {

        for (const auto &p : allPoints) {
            if (p.size() != allPoints[0].size()) {
                CurrentLogger::get()->error("Sparse points of differing dimensions", true);
            }
        }

        const vector<size_t> gridDimensions = _gridDimensions.value_or(
             vector<size_t>(
                    allPoints[0].size(),
                    ceil(pow(_nSparsePoints, 1.0 / static_cast<double>(allPoints[0].size())))
                    )
        );

        if (_minPoint.has_value() && _minPoint->size() != gridDimensions.size()) {
            CurrentLogger::get()->logAndThrow("Sparsifier's min_point dimensions {} don't match grid dimensions {}",
                                              _minPoint->size(), gridDimensions.size());
        }
        if (_maxPoint.has_value() && _maxPoint->size() != gridDimensions.size()) {
            CurrentLogger::get()->logAndThrow("Sparsifier's max_point dimensions {} don't match grid dimensions {}",
                                              _maxPoint->size(), gridDimensions.size());
        }

        vector<double> minPoint(gridDimensions.size());
        vector<double> maxPoint(gridDimensions.size());
        vector<double> step(gridDimensions.size());

        // TODO: looks ugly
        for (size_t d = 0; d < gridDimensions.size(); d++) {
            if (_minPoint.has_value()) {
                minPoint[d] = _minPoint.value()[d];
            } else {
                minPoint[d] = numeric_limits<double>::max();
                for (const auto &p : allPoints) {
                    minPoint[d] = min(minPoint[d], p[d]);
                }
            }
            if (_maxPoint.has_value()) {
                maxPoint[d] = _maxPoint.value()[d];
            } else {
                maxPoint[d] = numeric_limits<double>::min();
                for (const auto &p : allPoints) {
                    maxPoint[d] = max(maxPoint[d], p[d]+0.0001/*keep all points in bounds*/);
                }
            }
            step[d] = (maxPoint[d] - minPoint[d]) / static_cast<double>(gridDimensions[d]);
        }

        CurrentLogger::get()->info(
                "{}d histogram in range {} - {} with {} long bins:",
                gridDimensions.size(),
                iteratorToString(minPoint.begin(), minPoint.end()),
                iteratorToString(maxPoint.begin(), maxPoint.end()),
                iteratorToString(step.begin(), step.end())
                );

        vector<vector<double>> sparsePoints;
        map<vector<size_t>, vector<size_t>> usefulGridSlots;

        for (size_t i = 0; i < allPoints.size(); i++) {
            vector<size_t> gridSlot{};
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
        CurrentLogger::get()->debug(
            "Found {} grid slots containing some points -> attempting to sample each {} times / {} assigned randomly",
            usefulGridSlots.size(), reps, leftover
            );

        mt19937 gen(_seed);
        vector<vector<size_t>> usefulGridSlotsArr = {};
        for (auto &[gridSlot, pointIndices]: usefulGridSlots) {
            usefulGridSlotsArr.push_back(gridSlot);

            ranges::shuffle(pointIndices.begin(), pointIndices.end(), gen);
            for (size_t rep = 0; rep < min(pointIndices.size(), reps); rep++) {
                sparsePoints.push_back(allPoints[pointIndices[rep]]);
            }
        }

        uniform_int_distribution<> indexDist(0, usefulGridSlotsArr.size() - 1);

        while (sparsePoints.size() < _nSparsePoints) {
            const vector<size_t>& gridSlot = usefulGridSlotsArr[indexDist(gen)];

            vector<double> point(gridDimensions.size());
            for (size_t d = 0; d < gridDimensions.size(); d++) {
                uniform_real_distribution<> marginDist(0, step[d]);
                point[d] = minPoint[d] + step[d] * static_cast<double>(gridSlot[d]) + marginDist(gen);
            }

            sparsePoints.push_back(point);
            CurrentLogger::get()->debug(iteratorToString(point.begin(), point.end()));
        }

        vector<KernelParams> sparseKernelsParams;
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