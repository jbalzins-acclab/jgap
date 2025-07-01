#ifndef MATRIXBLOCK_HPP
#define MATRIXBLOCK_HPP

#include <utility>
#include <concepts>

#include "io/log/Logger.hpp"

using namespace std;

namespace jgap {

    class MatrixBlock {
    public:
        MatrixBlock(const size_t rows, const size_t columns)
            : _rows(rows), _columns(columns) {
            try {
                Logger::logger->info("Allocating matrix block");
                _data = vector<double>(rows * columns);
            } catch (const std::bad_alloc& e) {
                Logger::logger->error("Memory allocation failed");
            }
            Logger::logger->info("Matrix block allocated");
        };
        ~MatrixBlock() = default;

        double& operator()(const size_t i, const size_t j) {
            return _data[i * _columns + j];
        }

        vector<double>& rawData() {return _data;}

        [[nodiscard]] size_t rows() const {return _rows;}
        [[nodiscard]] size_t columns() const {return _columns;}
        [[nodiscard]] pair<size_t, size_t> blockSize() const {return{_rows, _columns};}

    private:
        size_t _rows;
        size_t _columns;
        vector<double> _data;
    };
}

#endif
