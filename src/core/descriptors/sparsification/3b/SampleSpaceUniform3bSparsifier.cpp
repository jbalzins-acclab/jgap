#include "core/descriptors/sparsification/3b/SampleSpaceUniform3bSparsifier.hpp"

#include <random>

#include "utils/Utils.hpp"

namespace jgap {
    SampleSpaceUniform3bSparsifier::SampleSpaceUniform3bSparsifier(const nlohmann::json &params)
        : PerSpecies3bSparsifier(params) {

        _nSparsePointsPerTriplet = params["n_sparse_per_triplet"].get<size_t>();

        if (params.contains("sparse_ranges")) {
            array<array<double, 2>, 3> sparseRanges{};
            sparseRanges[0][0] = params["sparse_ranges"]["+"]["from"].get<double>();
            sparseRanges[0][1] = params["sparse_ranges"]["+"]["to"].get<double>();
            sparseRanges[1][0] = params["sparse_ranges"]["--"]["from"].get<double>();
            sparseRanges[1][1] = params["sparse_ranges"]["--"]["to"].get<double>();
            sparseRanges[2][0] = params["sparse_ranges"]["jk"]["from"].get<double>();
            sparseRanges[2][1] = params["sparse_ranges"]["jk"]["to"].get<double>();
            _sparseRanges = sparseRanges;
        }

        if (params.contains("n_grid_points")) {
            _nGridPoints = array{
                params["n_grid_points"]["+"].get<size_t>(),
                params["n_grid_points"]["--"].get<size_t>(),
                params["n_grid_points"]["jk"].get<size_t>(),
            };
        }
    }

    PointsPerSpeciesTriplet SampleSpaceUniform3bSparsifier::sparsifyFromTripletsInData(PointsPerSpeciesTriplet& all3b) {
        PointsPerSpeciesTriplet result = {};

        for (auto &[speciesTriplet, tripletVectors] : all3b) {
            result[speciesTriplet] = {};

            Vector3 maxPoint = {0, 0, 0}, minPoint = {1e9, 1e9, 1e9}; // TODO: salkfjl

            for (auto &triplet : tripletVectors) {
                auto invariantTriplet = toInvariantTriplet({triplet.x, triplet.y}, triplet.z);

                minPoint.x = min(minPoint.x, invariantTriplet.x);
                minPoint.y = min(minPoint.y, invariantTriplet.y);
                minPoint.z = min(minPoint.z, invariantTriplet.z);

                maxPoint.x = max(maxPoint.x, invariantTriplet.x);
                maxPoint.y = max(maxPoint.y, invariantTriplet.y);
                maxPoint.z = max(maxPoint.z, invariantTriplet.z);
            }

            if (_sparseRanges.has_value()) {
                auto sr = _sparseRanges.value();
                minPoint = {sr[0][0], sr[1][0], sr[2][0]};
                maxPoint = {sr[0][1], sr[1][1], sr[2][1]};
            }

            const size_t n = _nSparsePointsPerTriplet;
            array<size_t, 3> nSteps = _nGridPoints.value_or(array{
                static_cast<size_t>(ceil(pow(n, 1.0/3.0))),
                static_cast<size_t>(ceil(pow(n, 1.0/3.0))),
                static_cast<size_t>(ceil(pow(n, 1.0/3.0)))
            });

            array steps = {
                Vector3{(maxPoint.x - minPoint.x) / static_cast<double>(nSteps[0]), 0, 0},
                Vector3{0, (maxPoint.y - minPoint.y) / static_cast<double>(nSteps[1]), 0},
                Vector3{0, 0, (maxPoint.z - minPoint.z) / static_cast<double>(nSteps[2])}
            };

            vector hist(nSteps[0], vector(nSteps[1], vector<size_t>(nSteps[2], 0)));
            vector<array<size_t, 3>> usefulIndexes;

            CurrentLogger::get()->info(format(
                "3b histogram of sides {},{},{} in range {} - {}",
                nSteps[0], nSteps[1], nSteps[2], minPoint.toString(), maxPoint.toString()
                ));
            CurrentLogger::get()->debug("3b " + speciesTriplet.toString() + " sparse points:");
            for (const Vector3 &originalPoint: tripletVectors) {
                auto point = toInvariantTriplet({originalPoint.x, originalPoint.y}, originalPoint.z);
                auto xIndex = static_cast<size_t>((point.x - minPoint.x) / steps[0].x);
                auto yIndex = static_cast<size_t>((point.y - minPoint.y) / steps[1].y);
                auto zIndex = static_cast<size_t>((point.z - minPoint.z) / steps[2].z);

                if (point.x >= maxPoint.x) xIndex--;
                if (point.y >= maxPoint.y) yIndex--;
                if (point.z >= maxPoint.z) zIndex--;

                if (++hist[xIndex][yIndex][zIndex] == 1) {
                    result[speciesTriplet].push_back(point);
                    CurrentLogger::get()->debug(point.toString());

                    usefulIndexes.push_back({xIndex, yIndex, zIndex});
                }
            }

            if (result[speciesTriplet].size() == n) {
                continue;
            }

            CurrentLogger::get()->debug("Not all 3b histogram bins have values => random selection:");

            mt19937 gen(799u); // NOLINT(*-msc51-cpp)
            uniform_real_distribution<> marginDistX(0, steps[0].x);
            uniform_real_distribution<> marginDistY(0, steps[1].y);
            uniform_real_distribution<> marginDistZ(0, steps[2].z);
            uniform_int_distribution<> indexDist(0, usefulIndexes.size() - 1);

            while (result[speciesTriplet].size() < n) {
                // select bin randomly
                const auto index = usefulIndexes[indexDist(gen)];

                result[speciesTriplet].push_back(
                    minPoint
                    + steps[0] * static_cast<double>(index[0])
                    + steps[1] * static_cast<double>(index[1])
                    + steps[2] * static_cast<double>(index[2])
                    + steps[0] * marginDistX(gen)
                    + steps[1] * marginDistY(gen)
                    + steps[2] * marginDistZ(gen)
                );
                CurrentLogger::get()->debug(result[speciesTriplet].back().toString());
            }
        }

        return result;
    }
}
