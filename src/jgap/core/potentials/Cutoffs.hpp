#ifndef JGAP_CUTOFFS_HPP
#define JGAP_CUTOFFS_HPP

#include <optional>
#include <map>
#include <algorithm>
#include <ranges>
#include <initializer_list>

#include "jgap/core/Real.hpp"
#include "jgap/utils/Utils.hpp"

namespace jgap {
    struct Cutoffs {
        std::map<size_t, Real> per_cluster_size;

        Cutoffs() = default;

        Cutoffs(std::map<size_t, Real> per_n_dependencies) : per_cluster_size(std::move(per_n_dependencies)) {}
        Cutoffs(const std::initializer_list<std::pair<const size_t, Real>> list) : per_cluster_size(list) {}
        Cutoffs(const std::pair<const size_t, Real> list) : per_cluster_size({list}) {}

        Real maxOverall() const {
            Real result = 0.0;
            for (const auto &cutoff: per_cluster_size | std::views::values) {
                result = std::max(result, static_cast<Real>(cutoff));
            }
            return result;
        }

        Real forDim(const size_t dim) const {
            if (per_cluster_size.contains(dim)) {
                return per_cluster_size.at(dim);
            }
            return 0.0;
        }

        Cutoffs operator+(const Cutoffs& other) const {
            auto res = *this;
            for (const auto& [deps, cutoff] : other.per_cluster_size) {
                if (res.per_cluster_size.contains(deps)) {
                    res.per_cluster_size[deps] = std::max(res.per_cluster_size[deps], cutoff);
                } else {
                    res.per_cluster_size[deps] = cutoff;
                }
            }
            return res;
        }

        Cutoffs& operator+=(const Cutoffs& other) {
            *this = *this + other;
            return *this;
        }
    };
}

#endif