#ifndef JGAP_UNSEQFOR_HPP
#define JGAP_UNSEQFOR_HPP

#include <algorithm>
#include <cstddef>
#include <utility>

#ifdef HAS_TBB
#include <memory>
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#endif

namespace jgap {

    /// @brief Limits the number of threads used by @ref unseqForEach and @ref unseqForIndex.
    class GlobalThreadingControl {
    public:
        explicit GlobalThreadingControl(size_t num_threads) {
#ifdef HAS_TBB
            if (num_threads > 0) {
                control_ =
                    std::make_unique<tbb::global_control>(tbb::global_control::max_allowed_parallelism, num_threads);
            }
#endif
        }

    private:
#ifdef HAS_TBB
        std::unique_ptr<tbb::global_control> control_;
#endif
    };

    /// @brief Apply func for each element in between [first, last) iterators, may be in parallel.
    /// @note Encapsulates tbb::parallel_for_each if TBB was linked successfully,
    /// but falls back to sequential for-each otherwise.
    /// This allows using parallelization inside the core without mandatory dependence on an external lib.
    template<typename Iterator, typename Func>
    void unseqForEach(Iterator first, Iterator last, Func&& func) {
#ifdef HAS_TBB
        tbb::parallel_for_each(first, last, std::forward<Func>(func));
#else
        std::for_each(first, last, std::forward<Func>(func));
#endif
    }

    /// @brief Apply func for each element of an iterable, may be in parallel.
    /// @note Encapsulates tbb::parallel_for_each if TBB was linked successfully,
    /// but falls back to sequential for-each otherwise.
    /// This allows using parallelization inside the core without mandatory dependence on an external lib.
    template<typename Iterable, typename Func>
    void unseqForEach(Iterable iterable, Func&& func) {
#ifdef HAS_TBB
        tbb::parallel_for_each(iterable.begin(), iterable.end(), std::forward<Func>(func));
#else
        std::for_each(iterable.begin(), iterable.end(), std::forward<Func>(func));
#endif
    }

    /// @brief Apply func for each index in range [first, last), may be in parallel.
    /// @note Encapsulates tbb::parallel_for if TBB was linked successfully,
    /// but falls back to sequential for-each otherwise.
    /// This allows using parallelization inside the core without mandatory dependence on an external lib.
    template<typename Func>
    void unseqForIndex(size_t first, size_t last, Func&& func) {
#ifdef HAS_TBB
        tbb::parallel_for(first, last, std::forward<Func>(func));
#else
        for (size_t i = first; i < last; ++i) {
            func(i);
        }
#endif
    }

}

#endif
