#ifndef JGAP_UNSEQFOR_HPP
#define JGAP_UNSEQFOR_HPP

#include <utility>
#include <algorithm>

#ifdef HAS_TBB
    #include <tbb/parallel_for_each.h>
    #include <tbb/parallel_for.h>
    #include <tbb/global_control.h>
    #include <memory>

    namespace jgap {
        
        class GlobalThreadingControl {
        public:
            explicit GlobalThreadingControl(size_t num_threads) {
                if (num_threads > 0) {
                    control_ = std::make_unique<tbb::global_control>(
                        tbb::global_control::max_allowed_parallelism, num_threads);
                }
            }
        private:
            std::unique_ptr<tbb::global_control> control_;
        };

        template <typename Iterator, typename Func>
        void unseqForEach(Iterator first, Iterator last, Func&& func) {
            tbb::parallel_for_each(first, last, std::forward<Func>(func));
        }

        template <typename Iterable, typename Func>
        void unseqForEach(Iterable iterable, Func&& func) {
            tbb::parallel_for_each(iterable.begin(), iterable.end(), std::forward<Func>(func));
        }

        template <typename Func>
        void unseqForIndex(size_t first, size_t last, Func&& func) {
            tbb::parallel_for(first, last, std::forward<Func>(func));
        }
    }
#else
    #include <cstddef>

    namespace jgap {

        class GlobalThreadingControl {
        public:
            explicit GlobalThreadingControl(size_t /*num_threads*/) {}
        };

        template <typename Iterator, typename Func>
        void unseqForEach(Iterator first, Iterator last, Func&& func) {
            std::for_each(first, last, std::forward<Func>(func));
        }

        template <typename Iterable, typename Func>
        void unseqForEach(Iterable iterable, Func&& func) {
            std::for_each(iterable.begin(), iterable.end(), std::forward<Func>(func));
        }

        template <typename Func>
        void unseqForIndex(size_t first, size_t last, Func&& func) {
            for (size_t i = first; i < last; ++i) {
                func(i);
            }
        }
    }
#endif

#endif
