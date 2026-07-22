#ifndef JGAP_FROMFILESPARSIFIER_HPP
#define JGAP_FROMFILESPARSIFIER_HPP

#include <fstream>
#include <vector>
#include "Sparsifier.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {

    template<size_t Dim>
    class FromFileSparsifier : public Sparsifier<Dim> {
    public:
        FromFileSparsifier(std::string filename) : filename(std::move(filename)) {}

        std::vector<Descriptor<Dim>> selectSparsePoints(
            const std::vector<Descriptor<Dim>>& descriptors) const override {
            std::ifstream file(filename);
            if (!file.is_open()) {
                JGAP_LOG_AND_THROW("Could not open file: {}", filename);
            }

            std::vector<Real> all_numbers;
            Real num;
            while (file >> num) {
                all_numbers.push_back(num);
            }

            if (all_numbers.size() % Dim != 0) {
                JGAP_LOG_AND_THROW(
                    "The number of values in the file ({}) is not divisible by the descriptor dimension ({})",
                    all_numbers.size(), Dim);
            }

            std::vector<Descriptor<Dim>> sparse_points(all_numbers.size() / Dim);
            for (size_t i = 0; i < all_numbers.size(); i += Dim) {
                for (size_t j = 0; j < Dim; ++j) {
                    sparse_points[i / Dim][j] = all_numbers[i * Dim + j];
                }
            }
            return sparse_points;
        }

    private:
        std::string filename;
    };
}

#endif
