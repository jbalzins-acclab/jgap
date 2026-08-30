#ifndef JGAP_QRKERNELFIT_HPP
#define JGAP_QRKERNELFIT_HPP

#include "../../../core/io/log/CurrentLogger.hpp"
#include "jgap/core/Matrix.hpp"
#include "jgap/core/fit/gap/GapFit.hpp"
#include "jgap/core/potentials/Potential.hpp"
#include "jgap/core/potentials/gap/GapPotential.hpp"

#include <memory>
#include <vector>

namespace jgap {

    class QRKernelFit : public GapFit {
    public:
        QRKernelFit() = default;

    protected:
        std::vector<Real> findCoefficients(
            std::vector<ValuePtr<GapComponent>>& gap_components,
            const std::vector<Atoms>& training_data,
            std::vector<EnergyData>& energies_without_external,
            std::vector<Regularization>& sigmas_inverse
        ) override;

        virtual std::vector<Real> leastSquares(Matrix<ColumnMajor>& A, std::vector<Real>& b);

    private:
        Matrix<ColumnMajor> formMatrixA(
            const std::vector<ValuePtr<GapComponent>>& gap_components,
            const std::vector<Atoms>& training_data,
            const std::vector<EnergyData>& energy_data,
            const std::vector<Regularization>& sigmas_inverse
        ) const;

        static std::vector<Real> formVectorB(
            const std::vector<ValuePtr<GapComponent>>& gap_components,
            const std::vector<EnergyData>& energy_data,
            const std::vector<Regularization>& sigmas_inverse
        );

        static void fillInverseSigmaLK_NM(
            const std::vector<ValuePtr<GapComponent>>& gap_components,
            const Atoms& atoms,
            const EnergyData& energy_data,
            const Regularization& sigmas_inverse,
            Matrix<ColumnMajor>& A,
            size_t starting_row
        );
    };
}

#endif
