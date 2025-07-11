#include "core/descriptors/sparsification/2b/SampleSpaceUniform2bSparsifier.hpp"

#include <random>

namespace jgap {
    SampleSpaceUniform2bSparsifier::SampleSpaceUniform2bSparsifier(const nlohmann::json &params)
        : PerSpecies2bSparsifier(params) {

        _nSparsePointsPerSpeciesPair = params["n_sparse_per_pair"];

        if (params.contains("sparse_range")) {
            _sparseRange = array{
                params["sparse_range"]["from"].get<double>(),
                params["sparse_range"]["to"].get<double>()
            };
        }
    }

    PointsPerSpeciesPair SampleSpaceUniform2bSparsifier::sparsifyFromDistancesInData(
                                                                PointsPerSpeciesPair& distancesInData) {

        PointsPerSpeciesPair result{};

        for (const auto& [speciesPair, distances] : distancesInData) {

            result[speciesPair] = vector<double>();

            double minDist, maxDist;
            if (_sparseRange.has_value()) {
                minDist = _sparseRange.value()[0];
                maxDist = _sparseRange.value()[1];
            } else {
                minDist = *ranges::min_element(distances);
                maxDist = *ranges::max_element(distances);
            }

            CurrentLogger::get()->info("2b sparse range = " + to_string(minDist) + "-" + to_string(maxDist));

            const auto n = _nSparsePointsPerSpeciesPair;

            const double step = (maxDist + 1e-6 - minDist) / static_cast<double>(n);
            vector<size_t> hist(n, 0);

            vector<size_t> usefulIndexes;

            CurrentLogger::get()->debug("2b " + speciesPair.toString() + " main sparse points:");
            for (double distance: distances) {
                const auto index = static_cast<size_t>((distance - minDist) / step);
                if (++hist[index] == 1) {
                    result[speciesPair].push_back(distance);
                    CurrentLogger::get()->debug(to_string(distance));
                    usefulIndexes.push_back(index);
                }
            }

            if (result[speciesPair].size() == n) {
                continue;
            }

            CurrentLogger::get()->debug("Not all 2b histogram bins have values => random selection:");

            mt19937 gen(9138741u);
            uniform_real_distribution<> marginDist(0, step);
            uniform_int_distribution<> indexDist(0, usefulIndexes.size() - 1);

            while (result[speciesPair].size() < n) {
                // select bin randomly
                const auto index = usefulIndexes[indexDist(gen)];
                result[speciesPair].push_back(
                    static_cast<double>(index) * step + marginDist(gen)
                );
                CurrentLogger::get()->debug(to_string(result[speciesPair].back()));
            }
        }

        return result;
    }
}
