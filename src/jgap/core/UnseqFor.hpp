#ifndef JGAP_UNSEQFOR_HPP
#define JGAP_UNSEQFOR_HPP

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <type_traits>
#include <utility>

#if defined(__APPLE__) || defined(HAS_GCD)
#include <dispatch/dispatch.h>
#elif defined(_OPENMP) || defined(HAS_OPENMP)
#include <omp.h>
#endif

namespace jgap {

    /// @brief Limits the number of threads used by @ref unseqForEach and @ref unseqForIndex.
    class GlobalThreadingControl {
    public:
        explicit GlobalThreadingControl(size_t num_threads) {
#if defined(__APPLE__) || defined(HAS_GCD)
            (void)num_threads;
#elif defined(_OPENMP) || defined(HAS_OPENMP)
            if (num_threads > 0) {
                omp_set_num_threads(static_cast<int>(num_threads));
            }
#else
            (void)num_threads;
#endif
        }
    };

    /// @brief Apply func for each element in between [first, last) iterators, may be in parallel.
    /// @note Uses Grand Central Dispatch on Apple platforms, OpenMP on platforms where OpenMP
    /// is enabled, and falls back to sequential std::for_each otherwise.
    template<typename Iterator, typename Func>
    void unseqForEach(Iterator first, Iterator last, Func&& func) {
#if defined(__APPLE__) || defined(HAS_GCD)
        if constexpr (
            std::is_base_of_v<
                std::random_access_iterator_tag,
                typename std::iterator_traits<Iterator>::iterator_category>
        ) {
            const auto n = last - first;
            if (n <= 0) {
                return;
            }
            if (n == 1) {
                func(*first);
                return;
            }
            struct Context {
                Iterator first_it;
                std::remove_reference_t<Func>* fn;
            } ctx{first, &func};

            dispatch_apply_f(
                static_cast<size_t>(n),
                DISPATCH_APPLY_AUTO,
                &ctx,
                [](void* context, size_t i) {
                    auto* c = static_cast<Context*>(context);
                    (*c->fn)(*(c->first_it + i));
                }
            );
        } else {
            std::for_each(first, last, std::forward<Func>(func));
        }
#elif defined(_OPENMP) || defined(HAS_OPENMP)
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
    /// @note Uses Grand Central Dispatch on Apple platforms, OpenMP on platforms where OpenMP
    /// is enabled, and falls back to sequential loop otherwise.
    template<typename Func>
    void unseqForIndex(size_t first, size_t last, Func&& func) {
        if (first >= last) {
            return;
        }
        const size_t n = last - first;
        if (n == 1) {
            func(first);
            return;
        }
#if defined(__APPLE__) || defined(HAS_GCD)
        struct Context {
            size_t start_idx;
            std::remove_reference_t<Func>* fn;
        } ctx{first, &func};

        dispatch_apply_f(
            n,
            DISPATCH_APPLY_AUTO,
            &ctx,
            [](void* context, size_t i) {
                auto* c = static_cast<Context*>(context);
                (*c->fn)(c->start_idx + i);
            }
        );
#elif defined(_OPENMP) || defined(HAS_OPENMP)
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
