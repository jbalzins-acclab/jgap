#ifndef LARGEMATRIX_HPP
#define LARGEMATRIX_HPP

#include <utility>
#include <concepts>
#include <stdexcept>

#include "memory/MatrixBlock.hpp"

using namespace std;

namespace jgap {

    class LargeMatrix {
    public:
        LargeMatrix(const pair<size_t, size_t> &nBlocks, const pair<size_t, size_t> &blockSize) {

            _blocks = vector(nBlocks.first, vector<shared_ptr<MatrixBlock>>(nBlocks.second));

            _blockRows = blockSize.first;
            _blockColumns = blockSize.second;
        }
        virtual ~LargeMatrix() = default;

        [[nodiscard]]
        pair<size_t, size_t> getNumberOfBlocks() const {
            return {_blocks.size(), _blocks[0].size()};
        }
        // virtual vector<pair<size_t, size_t>> blocks();

        virtual void loadBlock(size_t blockRow, size_t blockCol) = 0;

        shared_ptr<MatrixBlock> getBlock(size_t blockRow, size_t blockCol) {
            return _blocks[blockRow][blockCol];
        }

        void setBlock(const size_t blockRow, const size_t blockCol, shared_ptr<MatrixBlock> block) {
            _blocks[blockRow][blockCol] = std::move(block);
        }

        virtual void commit(size_t blockRow, size_t blockCol) {};

        void deloadBlock(size_t blockRow, size_t blockCol) {
            _blocks[blockRow][blockCol].reset();
        }

    protected:

        size_t _blockRows;
        size_t _blockColumns;

        vector<vector<shared_ptr<MatrixBlock>>> _blocks;
    };
}

#endif
