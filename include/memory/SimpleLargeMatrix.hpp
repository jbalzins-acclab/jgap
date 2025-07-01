//
// Created by Jegors Balzins on 18.6.2025.
//

#ifndef SIMPLELARGEMATRIX_HPP
#define SIMPLELARGEMATRIX_HPP

#include <utility>
#include <concepts>
#include <memory>
#include <stdexcept>
#include <vector>

#include "memory/LargeMatrix.hpp"
#include "memory/MatrixBlock.hpp"

using namespace std;

namespace jgap {

    class SimpleLargeMatrix : public LargeMatrix {
    public:
        explicit SimpleLargeMatrix(const shared_ptr<MatrixBlock>& block)
            : LargeMatrix({1, 1}, block->blockSize()) {
            _blocks[0][0] = block;
        }
        ~SimpleLargeMatrix() override = default;

        void loadBlock(size_t blockRow, size_t blockCol) override {};
        void commit(size_t blockRow, size_t blockCol) override {};
    };
}


#endif //SIMPLELARGEMATRIX_HPP
