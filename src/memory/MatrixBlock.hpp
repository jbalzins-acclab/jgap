#ifndef MATRIXBLOCK_HPP
#define MATRIXBLOCK_HPP

#include <vector>

#include "io/log/CurrentLogger.hpp"

namespace jgap {

    class MatrixBlock {
    public:
        MatrixBlock(const size_t rows, const size_t columns)
            : _rows(rows), _columns(columns) {
            try {
                _data = std::vector<double>(rows * columns);
            } catch (const std::bad_alloc&) {
                JGAP_LOG_AND_THROW("Block memory allocation failed");
            }
        }
        ~MatrixBlock() = default;

        double& operator()(const size_t i, const size_t j) {
            return _data[i * _columns + j];
        }

        std::vector<double>& rawData() {return _data;}

        [[nodiscard]] size_t rows() const {return _rows;}
        [[nodiscard]] size_t columns() const {return _columns;}
        [[nodiscard]] std::pair<size_t, size_t> blockSize() const {return{_rows, _columns};}

    private:
        size_t _rows;
        size_t _columns;
        std::vector<double> _data;
    };
}

#endif
