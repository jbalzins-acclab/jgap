#ifndef JGAP_QRGAPFIT_HPP
#define JGAP_QRGAPFIT_HPP

#include "jgap/core/Matrix.hpp"
#include "jgap/core/fit/gap/GapFit.hpp"
#include "jgap/core/potentials/Potential.hpp"
#include "jgap/core/potentials/gap/GapPotential.hpp"
#include "../../../core/io/log/CurrentLogger.hpp"

#include <memory>
#include <vector>

namespace jgap {

    class QRGapFit : public GapFit {
    public:
        explicit QRGapFit(const Real jitter = 1e-8) : jitter(jitter) {}

    protected:
        double jitter;

        std::vector<Real> findCoefficients(
            std::vector<ValuePtr<GapComponent>>& gap_components, const std::vector<Atoms>& training_data,
            std::vector<EnergyData>& energies_without_external, std::vector<Regularization>& sigmas_inverse
        ) override;

        virtual std::vector<Real> leastSquares(Matrix<ColumnMajor>& A, std::vector<Real>& b);

        Matrix<ColumnMajor> formMatrixA(
            const std::vector<ValuePtr<GapComponent>>& gap_components, const std::vector<Atoms>& training_data,
            const std::vector<EnergyData>& energy_data, const std::vector<Regularization>& sigmas_inverse
        ) const;

        static std::vector<Real> formVectorB(
            const std::vector<ValuePtr<GapComponent>>& gap_components, const std::vector<EnergyData>& energy_data,
            const std::vector<Regularization>& sigmas_inverse
        );

        void fillU_mm(
            size_t starting_row, size_t starting_col, const ValuePtr<GapComponent>& gap_component,
            Matrix<ColumnMajor>& A
        ) const;

        static void fillInverseSigmaK_nm(
            const std::vector<ValuePtr<GapComponent>>& gap_components, const Atoms& atoms,
            const EnergyData& energy_data, const Regularization& sigmas_inverse, Matrix<ColumnMajor>& A,
            size_t starting_row
        );

        static Matrix<ColumnMajor> choleskyDecomposition(Matrix<RowMajor>& matrix_block);
    };
}

#endif
