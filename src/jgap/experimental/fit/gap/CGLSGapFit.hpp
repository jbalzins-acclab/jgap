#ifndef JGAP_CGLSGAPFIT_HPP
#define JGAP_CGLSGAPFIT_HPP

#include <optional>
#include "jgap/core/Matrix.hpp"
#include "jgap/core/fit/gap/GapFit.hpp"

namespace jgap {

    class CGLSGapFit : public GapFit {
    public:
        explicit CGLSGapFit(
            const Real jitter = 1e-8, const std::optional<int> max_iterations = std::nullopt,
            const std::optional<double> tolerance = std::nullopt, const int print_interval = 1
        ) :
            jitter(jitter), max_iterations(max_iterations), tolerance(tolerance), print_interval(print_interval) {}

    protected:
        std::vector<Real> findCoefficients(
            std::vector<ValuePtr<GapComponent>>& gap_components, const std::vector<Atoms>& training_data,
            std::vector<EnergyData>& energies_without_external, std::vector<Regularization>& sigmas_inverse
        ) override;

        virtual std::vector<Real> leastSquares(Matrix<RowMajor>& A, std::vector<Real>& b);

    protected:
        Real jitter;
        std::optional<int> max_iterations;
        std::optional<double> tolerance;
        int print_interval{1};

        Matrix<RowMajor> formMatrixA(
            const std::vector<ValuePtr<GapComponent>>& gap_components, const std::vector<Atoms>& training_data,
            const std::vector<EnergyData>& energy_data, const std::vector<Regularization>& sigmas_inverse
        ) const;

        static std::vector<Real> formVectorB(
            const std::vector<ValuePtr<GapComponent>>& gap_components, const std::vector<EnergyData>& energy_data,
            const std::vector<Regularization>& sigmas_inverse
        );

        void fillU_mm(
            size_t starting_row, size_t starting_col, const ValuePtr<GapComponent>& gap_component, Matrix<RowMajor>& A
        ) const;

        static void fillInverseSigmaLK_NM(
            const std::vector<ValuePtr<GapComponent>>& gap_components, const Atoms& atoms,
            const EnergyData& energy_data, const Regularization& sigmas_inverse, Matrix<RowMajor>& A,
            size_t starting_row
        );

        static Matrix<RowMajor> choleskyDecomposition(Matrix<RowMajor>& matrix_block);
    };
}

#endif
