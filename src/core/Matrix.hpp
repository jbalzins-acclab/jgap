#ifndef MATRIXBLOCK_HPP
#define MATRIXBLOCK_HPP

#include <vector>

#include "Real.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {

    class Matrix {
    public:
        Matrix(const size_t rows, const size_t columns)
            : _rows(rows), _columns(columns) {
            try {
                _data = std::vector<Real>(rows * columns);
            } catch (const std::bad_alloc&) {
                JGAP_LOG_AND_THROW("Matrix memory allocation failed");
            }
        }
        ~Matrix() = default;

        Real& operator()(const size_t i, const size_t j) {
            return _data[i * _columns + j];
        }

        std::vector<Real>& rawData() {return _data;}

        size_t rows() const {return _rows;}
        size_t columns() const {return _columns;}
        std::pair<size_t, size_t> blockSize() const {return{_rows, _columns};}

    private:
        size_t _rows;
        size_t _columns;
        std::vector<Real> _data;
    };
}

#endif
