#ifndef JGAP_UNSEQFOR_HPP
#define JGAP_UNSEQFOR_HPP

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <type_traits>
#include <utility>

#if defined(_OPENMP) || defined(HAS_OPENMP)
#include <omp.h>
#endif

namespace jgap {

    /// @brief Limits the number of threads used by @ref unseqForEach and @ref unseqForIndex.
    class GlobalThreadingControl {
    public:
        explicit GlobalThreadingControl(size_t num_threads) {
#if defined(_OPENMP) || defined(HAS_OPENMP)
            if (num_threads > 0) {
                omp_set_num_threads(static_cast<int>(num_threads));
            }
#endif
        }
    };

    /// @brief Apply func for each element in between [first, last) iterators, may be in parallel.
    /// @note Uses OpenMP parallel for random-access iterators if OpenMP is enabled,
    /// falls back to sequential std::for_each otherwise.
    template<typename Iterator, typename Func>
    void unseqForEach(Iterator first, Iterator last, Func&& func) {
#if defined(_OPENMP) || defined(HAS_OPENMP)
        if constexpr (
            std::is_base_of_v<
                std::random_access_iterator_tag,
                typename std::iterator_traits<Iterator>::iterator_category>
        ) {
            const auto n = last - first;
            // clang-format off
            #pragma omp parallel for schedule(dynamic)
            // clang-format on
            for (auto i = 0; i < n; ++i) {
                func(*(first + i));
            }
        } else {
            std::for_each(first, last, std::forward<Func>(func));
        }
#else
        std::for_each(first, last, std::forward<Func>(func));
#endif
    }

    /// @brief Apply func for each element of an iterable, may be in parallel.
    template<typename Iterable, typename Func>
    void unseqForEach(Iterable&& iterable, Func&& func) {
        unseqForEach(std::begin(iterable), std::end(iterable), std::forward<Func>(func));
    }

    /// @brief Apply func for each index in range [first, last), may be in parallel.
    /// @note Uses OpenMP parallel loop if OpenMP is enabled,
    /// falls back to sequential loop otherwise.
    template<typename Func>
    void unseqForIndex(size_t first, size_t last, Func&& func) {
#if defined(_OPENMP) || defined(HAS_OPENMP)
        // clang-format off
        #pragma omp parallel for schedule(dynamic)
        // clang-format on
        for (size_t i = first; i < last; ++i) {
            func(i);
        }
#else
        for (size_t i = first; i < last; ++i) {
            func(i);
        }
#endif
    }

}

#endif
