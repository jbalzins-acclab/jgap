#ifndef JGAP_GRID1D_HPP
#define JGAP_GRID1D_HPP

#include <cassert>
#include <array>
#include <cstddef>
#include <iterator>
#include <utility>

#include "Vector3.hpp"
#include "io/log/CurrentLogger.hpp"

using namespace std;

namespace jgap {
    class Grid1d {
    public:
        Grid1d() = default;
        Grid1d(const Grid1d&) = default;  // allow copying
        Grid1d(Grid1d&&) noexcept = default;  // allow moving

        Grid1d(size_t size, double spacing, double origin)
            : spacing(spacing), origin(origin)
        {
            data.resize(size, 0.0);
        }

        Grid1d(vector<double> data, double spacing, double origin)
            : data(std::move(data)), spacing(spacing), origin(origin) {}

        // --- Element access ---
        double& operator()(size_t i) {
            assert(i < size());
            return data[i];
        }

        const double& operator()(size_t i) const {
            assert(i < size());
            return data[i];
        }

        Grid1d& operator=(const Grid1d&) = default;
        Grid1d& operator=(Grid1d&&) noexcept = default;

        Grid1d operator+(const Grid1d& other) const {
            checkShape(other);
            Grid1d result(size(), spacing, origin);
            for (size_t i = 0; i < size(); i++) {
                result.data[i] = data[i] + other.data[i];
            }
            return result;
        }

        Grid1d& operator+=(const Grid1d& other) {
            checkShape(other);
            for (size_t i = 0; i < size(); i++) {
                data[i] += other.data[i];
            }
            return *this;
        }

        // --- Coordinate conversion ---
        double coord(const size_t i) const {
            return origin + static_cast<double>(i) * spacing;
        }

        // --- Find closest <= grid index for a given coordinate ---
        size_t lowerIndex(double x) const {
            assert(x >= origin && "Point outside Grid1d(too low)");

            auto idx = static_cast<size_t>((x - origin) / spacing);
            assert(idx < size() && "Point outside Grid1d(too high)");

            return idx;
        }

        size_t size() const { return data.size(); }

        double spacing{};
        double origin{};
        vector<double> data{};

        void checkShape(const Grid1d& other) const {
            assert(size() == other.size() && "Grid1d dimensions don't match");
            assert(fabs(spacing - other.spacing) < 1e-12 && "Grid1d spacings don't match");
            assert(fabs(origin - other.origin) < 1e-12 && "Grid1d origins don't match");
        }

        // --- Iterator support ---
        struct CellRef {
            size_t index;
            double pos;
            double& value;
        };

        struct ConstCellRef {
            size_t index;
            double position;
            const double& value;
        };

        class Iterator {
        public:
            Iterator() : g(nullptr), i(0) {}
            Iterator(Grid1d* g, size_t i) : g(g), i(i) {}

            CellRef operator*() const {
                return CellRef{i, g->coord(i), g->data[i]};
            }

            Iterator& operator++() {
                ++i;
                return *this;
            }

            Iterator operator++(int) {
                Iterator tmp(*this);
                ++i;
                return tmp;
            }

            bool operator==(const Iterator& other) const {
                return g == other.g && i == other.i;
            }

            bool operator!=(const Iterator& other) const {
                return !(*this == other);
            }

        private:
            Grid1d* g;
            size_t i;
        };

        class ConstIterator {
        public:
            ConstIterator() : g(nullptr), i(0) {}
            ConstIterator(const Grid1d* g, size_t i) : g(g), i(i) {}

            ConstCellRef operator*() const {
                return ConstCellRef{i, g->coord(i), g->data[i]};
            }

            ConstIterator& operator++() {
                ++i;
                return *this;
            }

            ConstIterator operator++(int) {
                ConstIterator tmp(*this);
                ++i;
                return tmp;
            }

            bool operator==(const ConstIterator& other) const {
                return g == other.g && i == other.i;
            }

            bool operator!=(const ConstIterator& other) const {
                return !(*this == other);
            }

        private:
            const Grid1d* g;
            size_t i;
        };

        Iterator begin() { return Iterator(this, 0); }
        Iterator end() { return Iterator(this, size()); }
        ConstIterator begin() const { return ConstIterator(this, 0); }
        ConstIterator end() const { return ConstIterator(this, size()); }
        ConstIterator cbegin() const { return ConstIterator(this, 0); }
        ConstIterator cend() const { return ConstIterator(this, size()); }
    };

}

#endif