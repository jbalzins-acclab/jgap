#include "core/descriptors/sparsification/SampleSpaceUniformSparsifier.hpp"

namespace jgap {

    SampleSpaceUniformSparsifier::SampleSpaceUniformSparsifier(const nlohmann::json &params) {
        _nSparsePoints = params["n_sparse"].get<size_t>();
        _seed = params.value("seed", 799u);
        if (params.contains("grid_dimensions")) {
            _gridDimensions = vector<size_t>();
            for (const auto &dim : params["grid_dimensions"]) {
                _gridDimensions->push_back(dim.get<size_t>());
            }
        }
    }

    vector<vector<double>> SampleSpaceUniformSparsifier::selectSparsePoints(const vector<vector<double>> &allPoints) {

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

        vector<double> minPoint(gridDimensions.size());
        vector<double> maxPoint(gridDimensions.size());
        vector<double> step(gridDimensions.size());

        for (size_t d = 0; d < gridDimensions.size(); d++) {
            minPoint[d] = std::numeric_limits<double>::max();
            maxPoint[d] = 0;
            for (const auto &p : allPoints) {
                minPoint[d] = min(minPoint[d], p[d]);
                maxPoint[d] = max(maxPoint[d], p[d]+0.0001/*keep all points in bounds*/);
            }
            step[d] = (maxPoint[d] - minPoint[d]) / static_cast<double>(gridDimensions[d]);
        }

        CurrentLogger::get()->info(format(
                "{}d histogram in range {} - {} with {} long bins:",
                gridDimensions.size(),
                iteratorToString(minPoint.begin(), minPoint.end()),
                iteratorToString(maxPoint.begin(), maxPoint.end()),
                iteratorToString(step.begin(), step.end())
                ));

        vector<vector<double>> sparsePoints;
        set<vector<size_t>> usefulGridSlots;
        for (const auto &p : allPoints) {
            vector<size_t> gridSlot{};
            for (size_t d = 0; d < gridDimensions.size(); d++) {
                gridSlot.push_back((p[d] - minPoint[d]) / step[d]);
            }

            if (usefulGridSlots.insert(gridSlot).second) {
                sparsePoints.push_back(p);
                CurrentLogger::get()->debug(iteratorToString(p.begin(), p.end()));
            }
        }
        const vector usefulGridSlotsArr(usefulGridSlots.begin(), usefulGridSlots.end());

        if (sparsePoints.size() == _nSparsePoints) return sparsePoints;
        CurrentLogger::get()->debug("Not all histogram bins have values => random selection:");

        mt19937 gen(_seed);
        uniform_int_distribution<> indexDist(0, usefulGridSlotsArr.size() - 1);

        while (sparsePoints.size() < _nSparsePoints) {
            vector<size_t> gridSlot = usefulGridSlotsArr[indexDist(gen)];

            vector<double> point(gridDimensions.size());
            for (size_t d = 0; d < gridDimensions.size(); d++) {
                uniform_real_distribution<> marginDist(0, step[d]);
                point[d] = minPoint[d] + step[d] * static_cast<double>(gridSlot[d]) + marginDist(gen);
            }

            sparsePoints.push_back(point);
            CurrentLogger::get()->debug(iteratorToString(point.begin(), point.end()));
        }

        return sparsePoints;
    }
}