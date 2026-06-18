#ifndef JGAP_GRIDN_HPP
#define JGAP_GRIDN_HPP

#include <array>
#include <vector>
#include <cstddef>
#include <numeric>
#include <cassert>
#include <iterator>

#include "core/Real.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {

    template <size_t N>
    requires (N > 0)
    class Grid {
    public:
        std::array<size_t, N> dims{};
        std::array<Real, N> spacing{};
        std::array<Real, N> origin{};
        std::vector<Real> data_flat{};

        Grid() = default;
        Grid(const Grid&) = default;
        Grid(Grid&&) noexcept = default;
        Grid& operator=(const Grid&) = default;
        Grid& operator=(Grid&&) noexcept = default;

        Grid(const std::array<size_t, N>& dims,
             const std::array<Real, N>& spacing,
             const std::array<Real, N>& origin)
            : dims(dims), spacing(spacing), origin(origin) {

            size_t total_size = 1;
            for (size_t dim : dims) {
                total_size *= dim;
            }
            data_flat.resize(total_size);
            calculateStrides();
        }

        Grid(const std::array<size_t, N>& dims,
             const std::array<Real, N>& spacing,
             const std::array<Real, N>& origin,
             const std::vector<Real>& data_flat)
            : dims(dims), spacing(spacing), origin(origin), data_flat(data_flat) {

            size_t total_size = 1;
            for (size_t dim : dims) {
                total_size *= dim;
            }

            assert(data_flat.size() == total_size);

            calculateStrides();
        }

        Grid& operator+=(const Grid<N>& other) {
            check_shape(other);
            for (size_t i = 0; i < data_flat.size(); ++i) {
                data_flat[i] += other.data_flat[i];
            }
            return *this;
        }

        // --- Element access ---
        Real& operator()(const std::array<size_t, N>& indices) {
            return data_flat[getFlatIndex(indices)];
        }

        const Real& operator()(const std::array<size_t, N>& indices) const {
            return data_flat[getFlatIndex(indices)];
        }

        // --- Coordinate conversion ---
        std::array<Real, N> getCoord(const std::array<size_t, N>& indices) const {
            std::array<Real, N> coord;
            for (size_t i = 0; i < N; ++i) {
                coord[i] = origin[i] + static_cast<Real>(indices[i]) * spacing[i];
            }
            return coord;
        }

        std::array<size_t, N> getIndices(size_t flat_index) const {
            std::array<size_t, N> indices;
            for (size_t i = 0; i < N; ++i) {
                indices[i] = flat_index / strides[i];
                flat_index %= strides[i];
            }
            return indices;
        }

        std::array<Real, N> getCutoff() const {
            std::array<Real, N> res{};
            for (size_t i = 0; i < N; i++) {
                res[i] = origin[i] + static_cast<Real>(dims[i] - 1) * spacing[i];
            }
            return res;
        }

        // --- Find closest <= grid index for a given coordinate ---
        std::array<size_t, N> lowerIndex(const std::array<Real, N>& pos) const {
            std::array<size_t, N> indices;
            for (size_t i = 0; i < N; ++i) {
                assert(pos[i] >= origin[i] && "Point outside GridN (too low)");
                indices[i] = static_cast<size_t>((pos[i] - origin[i]) / spacing[i]);
                assert(indices[i] < dims[i] && "Point outside GridN (too high)");
            }
            return indices;
        }

        void checkShape(const Grid<N>& other) const {
            if (dims != other.dims) {
                JGAP_LOG_AND_THROW("GridN dimensions don't match");
            }
            for(size_t i = 0; i < N; ++i) {
                if (std::abs(spacing[i] - other.spacing[i]) > 1e-12) {
                    JGAP_LOG_AND_THROW("GridN spacings don't match");
                }
                if (std::abs(origin[i] - other.origin[i]) > 1e-12) {
                    JGAP_LOG_AND_THROW("GridN origins don't match");
                }
            }
        }

        // --- Iterator support ---
        struct CellRef {
            std::array<size_t, N> index;
            std::array<Real, N> pos;
            Real& value;
        };

        struct ConstCellRef {
            std::array<size_t, N> index;
            std::array<Real, N> pos;
            const Real& value;
        };

        template<bool IsConst>
        class base_iterator {
        public:
            using iterator_category = std::forward_iterator_tag;
            using grid_type = std::conditional_t<IsConst, const Grid<N>, Grid<N>>;
            using value_type = std::conditional_t<IsConst, ConstCellRef, CellRef>;
            using difference_type = std::ptrdiff_t;

            base_iterator() : g(nullptr), current_indices{} {}
            base_iterator(grid_type* g, std::array<size_t, N> indices) : g(g), current_indices(indices) {}

            value_type operator*() const {
                return {current_indices, g->getCoord(current_indices), (*g)(current_indices)};
            }

            base_iterator& operator++() {
                advance();
                return *this;
            }

            bool operator==(const base_iterator& other) const {
                return g == other.g && current_indices == other.current_indices;
            }
            bool operator!=(const base_iterator& other) const { return !(*this == other); }

        private:
            void advance() {
                if (!g) return;
                for (int i = N - 1; i >= 0; --i) {
                    if (++current_indices[i] < g->dims[i]) {
                        return;
                    }
                    current_indices[i] = 0;
                }
                // End of iteration, set to end state
                current_indices[0] = g->dims[0];
            }

            grid_type* g;
            std::array<size_t, N> current_indices;
        };

        using Iterator = base_iterator<false>;
        using ConstIterator = base_iterator<true>;

        Iterator begin() { return Iterator(this, {}); }
        Iterator end() {
            auto end_indices = std::array<size_t, N>{};
            if (dims[0] > 0) end_indices[0] = dims[0];
            return Iterator(this, end_indices);
        }
        ConstIterator begin() const { return ConstIterator(this, {}); }
        ConstIterator end() const {
            auto end_indices = std::array<size_t, N>{};
            if (dims[0] > 0) end_indices[0] = dims[0];
            return ConstIterator(this, end_indices);
        }
        ConstIterator cbegin() const { return begin(); }
        ConstIterator cend() const { return end(); }

    private:
        std::array<size_t, N> strides{};

        void calculateStrides() {
            strides[N - 1] = 1;
            for (size_t i = N - 1; i-- > 0;) {
                strides[i] = strides[i + 1] * dims[i + 1];
            }
        }

        size_t getFlatIndex(const std::array<size_t, N>& indices) const {
            size_t index = 0;
            for (size_t i = 0; i < N; ++i) {
                index += indices[i] * strides[i];
            }
            return index;
        }
    };
}

#endif