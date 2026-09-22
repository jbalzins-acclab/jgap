#ifndef JGAP_QRGAPFIT_HPP
#define JGAP_QRGAPFIT_HPP

#include "../../../core/io/log/CurrentLogger.hpp"
#include "jgap/core/Matrix.hpp"
#include "jgap/core/fit/gap/GapFit.hpp"
#include "jgap/core/potentials/Potential.hpp"
#include "jgap/core/potentials/gap/GapPotential.hpp"

#include <memory>
#include <vector>

namespace jgap {

    class QRGapFit : public GapFit {
    public:
        explicit QRGapFit(const double jitter = 1e-8) : jitter(jitter) {}

    protected:
        double jitter;

        void findCoefficients(
            GapPotential& to_be_fit,
            const std::vector<Atoms>& training_data,
            std::vector<EnergyData>& energies_without_external,
            std::vector<Regularization>& sigmas_inverse
        ) override;

        virtual std::vector<double> leastSquares(Matrix<ColumnMajor>& A, std::vector<double>& b);

        Matrix<ColumnMajor> formAugmentedCovarianceMatrixA(
            const std::vector<ValuePtr<GapComponent>>& gap_components,
            const std::vector<Atoms>& training_data,
            const std::vector<EnergyData>& energy_data,
            const std::vector<Regularization>& sigmas_inverse
        ) const;

        static std::vector<double> formTargetVectorBForEntry(
            const EnergyData& energy_data,
            const Regularization& sigmas_inverse
        );

        static std::vector<double> formNormalizedAugmentedTargetVectorB(
            const std::vector<ValuePtr<GapComponent>>& gap_components,
            const std::vector<EnergyData>& energy_data,
            const std::vector<Regularization>& sigmas_inverse
        );

        Matrix<ColumnMajor> U_MM(const ValuePtr<GapComponent>& gap_component) const;

        static Matrix<ColumnMajor> formInverseSigmaLK_NMForEntry(
            const std::vector<ValuePtr<GapComponent>>& gap_components,
            const Atoms& atoms,
            const EnergyData& energy_data,
            const Regularization& sigmas_inverse
        );

        static Matrix<ColumnMajor> choleskyDecomposition(Matrix<RowMajor>& matrix_block);
    };
}

#endif
