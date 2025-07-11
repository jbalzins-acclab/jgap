#include "core/descriptors/sparsification/eam/SampleSpaceUniformEamSparsifier.hpp"

#include <random>

using namespace std;

namespace jgap {
    SampleSpaceUniformEamSparsifier::SampleSpaceUniformEamSparsifier(const nlohmann::json& params)
            : PerSpeciesEamSparsifier(params) {

        _nSparsePoints = params["n_sparse_per_species"].get<size_t>();

        _sparseRanges = {};
        if (params.contains("sparse_ranges")) {
            for (const auto &[species, range]: params["sparse_ranges"].items()) {
                _sparseRanges[species] = {
                    range["from"].get<double>(),
                    range["to"].get<double>()
                };
            }
        }
    }

    DensitiesPerSpecies SampleSpaceUniformEamSparsifier::sparsifyFromDensitiesInData(
        DensitiesPerSpecies& densitiesInData) {

        DensitiesPerSpecies result{};

        for (const auto &[element, densities]: densitiesInData) {

            double rangeMin = 0, rangeMax = 0;
            if (!_sparseRanges.contains(element)) {
                rangeMin = ranges::min(densities);
                rangeMax = ranges::max(densities);
            } else {
                rangeMin = _sparseRanges[element][0];
                rangeMax = _sparseRanges[element][1];
            }

            const auto n = _nSparsePoints;

            CurrentLogger::get()->info(format(
                "EAM sparse range for {} = {}-{}",
                element,
                rangeMin,
                rangeMax
                ));
            const double step = (rangeMax - rangeMin) / static_cast<double>(n);

            vector<size_t> hist(n, 0);
            vector<size_t> usefulIndexes;

            CurrentLogger::get()->debug("EAM sparse points:");
            for (double density: densities) {
                const auto index = static_cast<size_t>((density - rangeMin) / step);

                //CurrentLogger::get()->debug(to_string(index) + "a");
                //CurrentLogger::get()->debug(to_string(n) + "b");
                if (++hist[index] == 1) {
                    result[element].push_back(density);
                    CurrentLogger::get()->debug(to_string(density));
                    usefulIndexes.push_back(index);
                }
            }

            if (result[element].size() == n) {
                continue;
            }

            CurrentLogger::get()->debug("Not all EAM histogram bins have values => random selection:");
            mt19937 gen(99);
            uniform_real_distribution<> marginDist(0, step);
            uniform_int_distribution<> indexDist(0, usefulIndexes.size() - 1);

            while (result[element].size() < n) {
                // select bin randomly
                const auto index = usefulIndexes[indexDist(gen)];
                result[element].push_back(
                    static_cast<double>(index) * step + marginDist(gen)
                );
                CurrentLogger::get()->debug(to_string(result[element].back()));
            }
        }
        return result;
    }
}
