#ifndef JGAP_HISTOGRAMUNIFORMSPARSIFIER_HPP
#define JGAP_HISTOGRAMUNIFORMSPARSIFIER_HPP

#include <cmath>
#include <map>
#include <random>
#include <set>

#include "Sparsifier.hpp"
#include "io/log/CurrentLogger.hpp"
#include "utils/Utils.hpp"

namespace jgap {

    template<size_t Dim>
    class HistogramUniformSparsifier : public Sparsifier<Dim> {
    public:
        HistogramUniformSparsifier(size_t seed, size_t n_sparse_points,
                                   std::optional<std::array<bool, Dim>> is_dim_active = std::nullopt,
                                   std::optional<std::array<size_t, Dim>> grid_dimensions = std::nullopt,
                                   std::optional<Descriptor<Dim>> min_point = std::nullopt,
                                   std::optional<Descriptor<Dim>> max_point = std::nullopt) :
            seed(seed),
            n_sparse_points(n_sparse_points),
            is_dim_active(std::move(is_dim_active)),
            grid_dimensions(std::move(grid_dimensions)),
            min_point(std::move(min_point)),
            max_point(std::move(max_point)) {}

        std::vector<Descriptor<Dim>> selectSparsePoints(const std::vector<Descriptor<Dim>>& descriptors) const override;

    private:
        size_t seed;
        size_t n_sparse_points;
        std::optional<std::array<bool, Dim>> is_dim_active;
        std::optional<std::array<size_t, Dim>> grid_dimensions;
        std::optional<Descriptor<Dim>> min_point;
        std::optional<Descriptor<Dim>> max_point;
    };

    template<size_t Dim>
    inline std::vector<Descriptor<Dim>> HistogramUniformSparsifier<Dim>::selectSparsePoints(
        const std::vector<Descriptor<Dim>>& descriptors) const {
        if (descriptors.empty()) {
            return {};
        }

        std::array<size_t, Dim> grid_dimensions_{};
        if (grid_dimensions.has_value()) {
            grid_dimensions_ = grid_dimensions.value();
        } else {
            size_t active_dims = 0;
            if (is_dim_active.has_value()) {
                assert(is_dim_active.value().size() == dim);
                for (bool active: is_dim_active.value()) {
                    if (active) active_dims++;
                }
            } else {
                active_dims = dim;
            }

            size_t active_grid_dim =
                active_dims > 0
                    ? static_cast<size_t>(std::ceil(std::pow(n_sparse_points, 1.0 / static_cast<Real>(active_dims))))
                    : 1;

            for (size_t d = 0; d < Dim; d++) {
                bool active = is_dim_active.has_value() ? is_dim_active.value()[d] : true;
                grid_dimensions_[d] = active ? active_grid_dim : 1;
            }
        }

        std::array<Real, Dim> min_point_{};
        std::array<Real, Dim> max_point_{};
        std::array<Real, Dim> step{};

        for (size_t d = 0; d < Dim; d++) {
            if (min_point.has_value()) {
                min_point_[d] = min_point.value()[d];
            } else {
                min_point_[d] = std::numeric_limits<Real>::max();
                for (size_t i = 0; i < descriptors.size(); ++i) {
                    min_point_[d] = std::min(min_point_[d], descriptors[i][d]);
                }
            }
            if (max_point.has_value()) {
                max_point_[d] = max_point.value()[d];
            } else {
                max_point_[d] = std::numeric_limits<Real>::lowest();
                for (size_t i = 0; i < descriptors.size(); ++i) {
                    max_point_[d] = std::max(max_point_[d], descriptors[i][d] + 0.0001 /*keep all points in bounds*/);
                }
            }
            step[d] = (max_point_[d] - min_point_[d]) / static_cast<Real>(grid_dimensions_[d]);
        }

        JGAP_LOG_INFO("{}d histogram in range {} - {} with {} long bins:", Dim,
                      iteratorToString(min_point_.begin(), min_point_.end()),
                      iteratorToString(max_point_.begin(), max_point_.end()),
                      iteratorToString(step.begin(), step.end()));

        std::vector<Descriptor<Dim>> sparse_points;
        std::map<std::array<size_t, Dim>, std::vector<size_t>> useful_grid_slots;

        for (size_t i = 0; i < descriptors.size(); i++) {
            std::array<size_t, Dim> grid_slot{};
            for (size_t d = 0; d < Dim; d++) {
                grid_slot[d] = static_cast<size_t>((descriptors[i][d] - min_point_[d]) / step[d]);
            }

            if (!useful_grid_slots.contains(grid_slot)) {
                useful_grid_slots[grid_slot] = {};
            }
            useful_grid_slots[grid_slot].push_back(i);
        }

        const size_t reps = n_sparse_points / useful_grid_slots.size();
        const size_t leftover = n_sparse_points - useful_grid_slots.size() * reps;

        JGAP_LOG_DEBUG(
            "Found {} grid slots containing some points -> attempting to sample each {} times / {} assigned "
            "randomly",
            useful_grid_slots.size(), reps, leftover);

        std::mt19937 gen(seed);
        std::vector<std::array<size_t, Dim>> useful_grid_slots_arr = {};
        for (auto& [grid_slot, point_indices]: useful_grid_slots) {
            useful_grid_slots_arr.push_back(grid_slot);

            std::ranges::shuffle(point_indices.begin(), point_indices.end(), gen);
            for (size_t rep = 0; rep < std::min(point_indices.size(), reps); rep++) {
                sparse_points.push_back(descriptors[point_indices[rep]]);
            }
        }

        std::uniform_int_distribution<> index_dist(0, (int) useful_grid_slots_arr.size() - 1);

        while (sparse_points.size() < n_sparse_points) {
            const std::array<size_t, Dim>& grid_slot = useful_grid_slots_arr[index_dist(gen)];

            Descriptor<Dim> pt;
            for (size_t d = 0; d < Dim; d++) {
                std::uniform_real_distribution<> margin_dist(0, step[d]);
                pt[d] = min_point_[d] + step[d] * static_cast<Real>(grid_slot[d]) + margin_dist(gen);
            }
            sparse_points.push_back(pt);
        }

        return sparse_points;
    }
}

#endif
