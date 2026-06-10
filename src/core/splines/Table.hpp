#ifndef JGAP_GRIDN_HPP
#define JGAP_GRIDN_HPP

#include <array>
#include <vector>
#include <cstddef>
#include <numeric>
#include <cassert>
#include <iterator>

#include "core/Real.hpp"

namespace jgap {

template <size_t N>
requires (N > 0)
class Table {
public:
    std::array<size_t, N> dims{};
    std::array<Real, N> spacing{};
    std::array<Real, N> origin{};
    std::vector<Real> data_flat{};

    Table() = default;
    Table(const Table&) = default;
    Table(Table&&) noexcept = default;
    Table& operator=(const Table&) = default;
    Table& operator=(Table&&) noexcept = default;

    Table(const std::array<size_t, N>& dims, const std::array<Real, N>& spacing, const std::array<Real, N>& origin)
        : dims(dims), spacing(spacing), origin(origin) {
        size_t total_size = 1;
        for (size_t dim : dims) {
            total_size *= dim;
        }
        data_flat.resize(total_size);
        calculate_strides();
    }

    // --- Element access ---
    Real& operator()(const std::array<size_t, N>& indices) {
        return data_flat[get_flat_index(indices)];
    }

    const Real& operator()(const std::array<size_t, N>& indices) const {
        return data_flat[get_flat_index(indices)];
    }

    // --- Coordinate conversion ---
    std::array<Real, N> get_coord(const std::array<size_t, N>& indices) const {
        std::array<Real, N> coord;
        for (size_t i = 0; i < N; ++i) {
            coord[i] = origin[i] + static_cast<Real>(indices[i]) * spacing[i];
        }
        return coord;
    }

    // --- Find closest <= grid index for a given coordinate ---
    std::array<size_t, N> lower_index(const std::array<Real, N>& pos) const {
        std::array<size_t, N> indices;
        for (size_t i = 0; i < N; ++i) {
            assert(pos[i] >= origin[i] && "Point outside GridN (too low)");
            indices[i] = static_cast<size_t>((pos[i] - origin[i]) / spacing[i]);
            assert(indices[i] < dims[i] && "Point outside GridN (too high)");
        }
        return indices;
    }

    void check_shape(const Table<N>& other) const {
        assert(dims == other.dims && "GridN dimensions don't match");
        for(size_t i = 0; i < N; ++i) {
            assert(std::abs(spacing[i] - other.spacing[i]) < 1e-12 && "GridN spacings don't match");
            assert(std::abs(origin[i] - other.origin[i]) < 1e-12 && "GridN origins don't match");
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
        using grid_type = std::conditional_t<IsConst, const Table<N>, Table<N>>;
        using value_type = std::conditional_t<IsConst, ConstCellRef, CellRef>;
        using difference_type = std::ptrdiff_t;

        base_iterator() : g(nullptr), current_indices{} {}
        base_iterator(grid_type* g, std::array<size_t, N> indices) : g(g), current_indices(indices) {}

        value_type operator*() const {
            return {current_indices, g->get_coord(current_indices), (*g)(current_indices)};
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

    void calculate_strides() {
        strides[N - 1] = 1;
        for (size_t i = N - 1; i-- > 0;) {
            strides[i] = strides[i + 1] * dims[i + 1];
        }
    }

    size_t get_flat_index(const std::array<size_t, N>& indices) const {
        size_t index = 0;
        for (size_t i = 0; i < N; ++i) {
            index += indices[i] * strides[i];
        }
        return index;
    }
};

}

#endif