#ifndef HISTOGRAMUNIFORMSPARSIFIER_HPP
#define HISTOGRAMUNIFORMSPARSIFIER_HPP

#include <cmath>
#include <random>
#include <set>
#include <numeric> // For std::iota

#include "Sparsifier.hpp"
#include "utils/Utils.hpp"
#include "io/log/CurrentLogger.hpp"
#include "../descriptors/Descriptor.hpp"

namespace jgap {

    template<size_t N_DIMENSIONS, size_t N_ATOMS>
    class HistogramUniformSparsifier : public Sparsifier<N_DIMENSIONS, N_ATOMS> {
    public:
        //SETUP_PARSER(Sparsifier, HistogramUniformSparsifier, histogram_uniform)

        //using Sparsifier<N_DIMENSIONS>::DescriptorValue;

        HistogramUniformSparsifier(size_t seed,
                                   size_t n_sparse_points,
                                   std::optional<std::array<size_t, N_DIMENSIONS>> grid_dimensions = std::nullopt,
                                   std::optional<std::array<double, N_DIMENSIONS>> min_point = std::nullopt,
                                   std::optional<std::array<double, N_DIMENSIONS>> max_point = std::nullopt)
            : seed_(seed),
              n_sparse_points_(n_sparse_points),
              grid_dimensions_(grid_dimensions),
              min_point_(min_point),
              max_point_(max_point) {
        }

        std::vector<Descriptor<N_DIMENSIONS, N_ATOMS>> selectSparsePoints(
            const std::vector<Descriptor<N_DIMENSIONS, N_ATOMS>> &descriptors);

        //std::vector<DescriptorValue<N_DIMENSIONS>> selectSparsePoints(const std::vector<std::vector<double>> &allPoints) override;

    private:
        size_t seed_;
        size_t n_sparse_points_;
        std::optional<std::array<size_t, N_DIMENSIONS>> grid_dimensions_;
        std::optional<std::array<double, N_DIMENSIONS>> min_point_;
        std::optional<std::array<double, N_DIMENSIONS>> max_point_;
    };

    template<size_t N_DIMENSIONS, size_t N_ATOMS>
    std::vector<Descriptor<N_DIMENSIONS, N_ATOMS>> HistogramUniformSparsifier<N_DIMENSIONS, N_ATOMS>
        ::selectSparsePoints(const std::vector<Descriptor<N_DIMENSIONS, N_ATOMS>> &descriptors) {

        if (descriptors.empty()) {
            return {};
        }

        std::array<size_t, N_DIMENSIONS> grid_dimensions;
        if (grid_dimensions_.has_value()) {
            grid_dimensions = grid_dimensions_.value();
        } else {
            grid_dimensions.fill(ceil(pow(n_sparse_points_, 1.0 / static_cast<double>(N_DIMENSIONS))));
        }

        std::array<double, N_DIMENSIONS> min_point{};
        std::array<double, N_DIMENSIONS> max_point{};
        std::array<double, N_DIMENSIONS> step{};

        for (size_t d = 0; d < N_DIMENSIONS; d++) {
            if (min_point_.has_value()) {
                min_point[d] = min_point_.value()[d];
            } else {
                min_point[d] = std::numeric_limits<double>::max();
                for (const auto &p : descriptors) {
                    min_point[d] = std::min(min_point[d], p.value[d]);
                }
            }
            if (max_point_.has_value()) {
                max_point[d] = max_point_.value()[d];
            } else {
                max_point[d] = std::numeric_limits<double>::min();
                for (const auto &p : descriptors) {
                    max_point[d] = std::max(max_point[d], p.value[d] + 0.0001/*keep all points in bounds*/);
                }
            }
            step[d] = (max_point[d] - min_point[d]) / static_cast<double>(grid_dimensions[d]);
        }

        JGAP_LOG_INFO(
            "{}d histogram in range {} - {} with {} long bins:",
            N_DIMENSIONS,
            iteratorToString(min_point.begin(), min_point.end()),
            iteratorToString(max_point.begin(), max_point.end()),
            iteratorToString(step.begin(), step.end())
        );

        std::vector<Descriptor<N_DIMENSIONS, N_ATOMS>> sparse_points;
        std::map<std::array<size_t, N_DIMENSIONS>, std::vector<size_t>> useful_grid_slots;

        for (size_t i = 0; i < descriptors.size(); i++) {
            std::array<size_t, N_DIMENSIONS> grid_slot{};
            for (size_t d = 0; d < N_DIMENSIONS; d++) {
                grid_slot[d] = static_cast<size_t>((descriptors[i].value[d] - min_point[d]) / step[d]);
            }

            if (!useful_grid_slots.contains(grid_slot)) {
                useful_grid_slots[grid_slot] = {};
            }
            useful_grid_slots[grid_slot].push_back(i);
        }
        const size_t reps = n_sparse_points_ / useful_grid_slots.size();
        const size_t leftover = n_sparse_points_ - useful_grid_slots.size() * reps;
        JGAP_LOG_DEBUG(
            "Found {} grid slots containing some points -> attempting to sample each {} times / {} assigned randomly",
            useful_grid_slots.size(), reps, leftover
        );

        std::mt19937 gen(seed_);
        std::vector<std::array<size_t, N_DIMENSIONS>> useful_grid_slots_arr = {};
        for (auto &[grid_slot, point_indices]: useful_grid_slots) {
            useful_grid_slots_arr.push_back(grid_slot);

            std::ranges::shuffle(point_indices.begin(), point_indices.end(), gen);
            for (size_t rep = 0; rep < std::min(point_indices.size(), reps); rep++) {
                sparse_points.push_back(descriptors[point_indices[rep]]);
            }
        }

        std::uniform_int_distribution<> index_dist(0, useful_grid_slots_arr.size() - 1);

        while (sparse_points.size() < n_sparse_points_) {
            const std::array<size_t, N_DIMENSIONS>& grid_slot = useful_grid_slots_arr[index_dist(gen)];

            Descriptor<N_DIMENSIONS, N_ATOMS> d{};
            for (size_t dim = 0; dim < N_DIMENSIONS; dim++) {
                std::uniform_real_distribution<> margin_dist(0, step[dim]);
                d.value[dim] = min_point[dim] + step[dim] * static_cast<double>(grid_slot[dim]) + margin_dist(gen);
            }

            sparse_points.push_back(d);
            JGAP_LOG_DEBUG(iteratorToString(d.value.begin(), d.value.end()));
        }

        return sparse_points;
    }
}

#endif
