#ifndef MATRIXBLOCK_HPP
#define MATRIXBLOCK_HPP

#include <utility>
#include <concepts>

#include "io/log/CurrentLogger.hpp"

using namespace std;

namespace jgap {

    class MatrixBlock {
    public:
        MatrixBlock(const size_t rows, const size_t columns)
            : _rows(rows), _columns(columns) {
            try {
                CurrentLogger::get()->debug(format ("Allocating {}x{} matrix block", rows, columns));
                _data = vector<double>(rows * columns);
            } catch (const std::bad_alloc& e) {
                CurrentLogger::get()->error("Memory allocation failed");
            }
            CurrentLogger::get()->info("Matrix block allocated");
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
