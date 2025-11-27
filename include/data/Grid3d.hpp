#ifndef JGAP_GRID3D_HPP
#define JGAP_GRID3D_HPP

#include <cassert>
#include <array>
#include <cstddef>
#include <iterator>

#include "Vector3.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {

    class Grid3d {
    public:
        size_t nR{}, nAngular{};
        Vector3 spacing{};
        Vector3 origin{};
        vector<double> dataFlat{};

        Grid3d() = default;
        Grid3d(const Grid3d&) = default;  // allow copying
        Grid3d(Grid3d&&) noexcept = default;  // allow moving
        Grid3d& operator=(const Grid3d&) = default;
        Grid3d& operator=(Grid3d&&) noexcept = default;

        Grid3d(double rMin, size_t nR, double spacingR, double angularMin, size_t nAngular, double spacingAngular)
            : nR(nR), nAngular(nAngular),
              spacing({spacingR, spacingR, spacingAngular}),
              origin({rMin, rMin, angularMin})
        {
            dataFlat.resize(nR * nR * nAngular);
        }

        double spacingR() const { return spacing.x; }
        double spacingAngular() const { return spacing.z; }
        double originR() const { return origin.x; }
        double originAngular() const { return origin.z; }

        double cutoff() const { return originR() + spacingR() * static_cast<double>(nR); }

        // --- Element access ---
        double& operator()(size_t i, size_t j, size_t k) {
            return dataFlat[i*nR*nR + j*nR + k];
        }

        const double& operator()(size_t i, size_t j, size_t k) const {
            return dataFlat[i*nR*nR + j*nR + k];
        }

        Grid3d operator+(const Grid3d& other) const {
            checkShape(other);
            Grid3d result(origin.x, nR, spacing.x, origin.z, nAngular, spacing.z);
            for (size_t i = 0; i < nR; i++) {
                for (size_t j = 0; j < nR; j++) {
                    for (size_t k = 0; k < nAngular; k++) {
                        result(i, j, k) = (*this)(i, j, k) + other(i, j, k);
                    }
                }
            }
            return result;
        }

        Grid3d& operator+=(const Grid3d& other) {
            checkShape(other);
            for (size_t i = 0; i < nR; i++) {
                for (size_t j = 0; j < nR; j++) {
                    for (size_t k = 0; k < nAngular; k++) {
                        (*this)(i, j, k) += other(i, j, k);
                    }
                }
            }
            return *this;
        }

        vector<double> slice_ij(size_t i, size_t j) const {
            vector<double> result(nAngular);
            for (size_t k = 0; k < nAngular; k++) {
                result[k] = dataFlat[i*nR*nR + j*nR + k];
            }
            return result;
        }

        vector<double> slice_ik(size_t i, size_t k) const {
            vector<double> result(nR);
            for (size_t j = 0; j < nR; j++) {
                result[j] = dataFlat[i*nR*nR + j*nR + k];
            }
            return result;
        }

        vector<double> slice_jk(size_t j, size_t k) const {
            vector<double> result(nR);
            for (size_t i = 0; i < nR; i++) {
                result[j] = dataFlat[i*nR*nR + j*nR + k];
            }
            return result;
        }

        // --- Coordinate conversion ---
        Vector3 coord(size_t i, size_t j, size_t k) const {
            return {
                origin.x + static_cast<double>(i) * spacing.x,
                origin.y + static_cast<double>(j) * spacing.y,
                origin.z + static_cast<double>(k) * spacing.z
            };
        }

        // --- Find closest <= grid index for a given coordinate ---
        array<size_t, 3> lowerIndex(const Vector3& pos) const {

            assert(pos.x >= origin.x && pos.y >= origin.y && pos.z >= origin.z && "Point outside Grid3d(too low)");

            auto idx = array{
                static_cast<size_t>((pos.x - origin.x) / spacing.x),
                static_cast<size_t>((pos.y - origin.y) / spacing.y),
                static_cast<size_t>((pos.z - origin.z) / spacing.z),
            };

            assert(idx[0] < nR && idx[1] < nR && idx[2] < nAngular && "Point outside Grid3d(too high)");

            return idx;
        }

        void checkShape(const Grid3d& other) const {
            assert(nR != other.nR || nR != other.nR || nAngular != other.nAngular && "Grid3d dimensions don't match");
            assert((spacing - other.spacing).len() < 1e-12 && "Grid3d spacings don't match");
            assert((origin - other.origin).len() < 1e-12 && "Grid3d origins don't match");
        }

        // --- Iterator support ---
        struct CellRef {
            array<size_t, 3> index;
            Vector3 pos;
            double& value;
        };
        struct ConstCellRef {
            array<size_t, 3> index;
            Vector3 position;
            const double& value;
        };

        class Iterator {
        public:
            using iterator_category = std::forward_iterator_tag;
            using value_type = CellRef;
            using difference_type = std::ptrdiff_t;
            using pointer = value_type*;
            using reference = value_type&;

            Iterator() : g(nullptr), i(0), j(0), k(0) {}
            Iterator(Grid3d* g, size_t i, size_t j, size_t k) : g(g), i(i), j(j), k(k) {}

            CellRef operator*() const {
                return CellRef{{i, j, k}, g->coord(i, j, k), (*g)(i, j ,k)};
            }

            Iterator& operator++() {
                advance();
                return *this;
            }
            Iterator operator++(int) {
                Iterator tmp(*this);
                advance();
                return tmp;
            }

            bool operator==(const Iterator& other) const {
                return g == other.g && i == other.i && j == other.j && k == other.k;
            }
            bool operator!=(const Iterator& other) const { return !(*this == other); }

        private:
            void advance() {
                if (!g) return;
                if (++k < g->nAngular) return;
                k = 0;
                if (++j < g->nR) return;
                j = 0;
                if (++i < g->nR) return;
                // end condition
                i = g->nR; j = 0; k = 0;
            }

            Grid3d* g;
            size_t i, j, k;
        };

        class ConstIterator {
        public:
            using iterator_category = std::forward_iterator_tag;
            using value_type = ConstCellRef;
            using difference_type = std::ptrdiff_t;
            using pointer = value_type*;
            using reference = value_type&;

            ConstIterator() : g(nullptr), i(0), j(0), k(0) {}
            ConstIterator(const Grid3d* g, size_t i, size_t j, size_t k) : g(g), i(i), j(j), k(k) {}

            ConstCellRef operator*() const {
                return ConstCellRef{{i, j, k}, g->coord(i, j, k), (*g)(i, j, k)};
            }

            ConstIterator& operator++() {
                advance();
                return *this;
            }
            ConstIterator operator++(int) {
                ConstIterator tmp(*this);
                advance();
                return tmp;
            }

            bool operator==(const ConstIterator& other) const {
                return g == other.g && i == other.i && j == other.j && k == other.k;
            }
            bool operator!=(const ConstIterator& other) const { return !(*this == other); }

        private:
            void advance() {
                if (!g) return;
                if (++k < g->nAngular) return;
                k = 0;
                if (++j < g->nR) return;
                j = 0;
                if (++i < g->nR) return;
                // end condition
                i = g->nR; j = 0; k = 0;
            }

            const Grid3d* g;
            size_t i, j, k;
        };

        Iterator begin() { return Iterator(this, 0, 0, 0); }
        Iterator end() { return Iterator(this, nR, 0, 0); }
        ConstIterator begin() const { return ConstIterator(this, 0, 0, 0); }
        ConstIterator end() const { return ConstIterator(this, nR, 0, 0); }
        ConstIterator cbegin() const { return ConstIterator(this, 0, 0, 0); }
        ConstIterator cend() const { return ConstIterator(this, nR, 0, 0); }
    };
}

#endif