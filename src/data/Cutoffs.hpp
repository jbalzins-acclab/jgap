#ifndef JGAP_CUTOFFS_HPP
#define JGAP_CUTOFFS_HPP

#include <optional>
#include <map>
#include <algorithm>
#include <ranges>

#include "core/Real.hpp"
#include "utils/Utils.hpp"

namespace jgap {
    struct Cutoffs {
        std::map<size_t, Real> per_n_dependencies;

        Cutoffs() = default;

        Cutoffs(std::map<size_t, Real> per_n_dependencies) : per_n_dependencies(std::move(per_n_dependencies)) {}

        Real maxOverall() const {
            Real result = 0.0;
            for (const auto &cutoff: per_n_dependencies | std::views::values) {
                result = std::max(result, static_cast<Real>(cutoff));
            }
            return result;
        }

        Cutoffs operator+(const Cutoffs& other) const {
            auto res = *this;
            for (const auto& [deps, cutoff] : other.per_n_dependencies) {
                if (res.per_n_dependencies.contains(deps)) {
                    res.per_n_dependencies[deps] = std::max(res.per_n_dependencies[deps], cutoff);
                } else {
                    res.per_n_dependencies[deps] = cutoff;
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