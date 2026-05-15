#ifndef JGAP_QRGAPFIT_HPP
#define JGAP_QRGAPFIT_HPP

#include "io/log/CurrentLogger.hpp"
#include "core/potentials/Potential.hpp"
#include "../../potentials/gap/GapPotential.hpp"

#include <Eigen/Dense>

#include <memory>

#include "GapFit.hpp"

namespace jgap {

    using EigenMatrix = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;
    using EigenVector = Eigen::Matrix<Real, Eigen::Dynamic, 1>;

    class QRGapFit : public GapFit {
    public:
        QRGapFit(Real jitter = 1e-8);

    protected:
        std::vector<Real> mainFit(
            std::vector<IGapComponent::Ptr> &gap_components,
            std::vector<NeighbourList> &neighbour_lists,
            std::vector<EnergyData> energies_without_external,
            std::vector<Regularization> &sigmas_inverse) override;

        std::vector<Real> leastSquares(EigenMatrix &A, EigenVector &b);

    private:
        double jitter;

        // TODO: 4-element matrix has only ~12% non-zero values i.e. ~12GB
        // => try to configure sparse matrix storage + linalg
        // * either specify manually to use it, or auto-select if #elements > 2 (?)

        EigenMatrix formMatrixA(const std::vector<IGapComponent::Ptr> &gap_components,
                                const std::vector<NeighbourList> &neighbour_lists,
                                const std::vector<EnergyData> &energy_data,
                                const std::vector<Regularization> &sigmas_inverse) const;

        static EigenVector formVectorB(const std::vector<IGapComponent::Ptr> &gap_components,
                                       const std::vector<EnergyData> &energy_data,
                                       const std::vector<Regularization> &sigmas_inverse);

        void fillU_mm(const size_t startingRow,
                      const size_t startingCol,
                      const IGapComponent::Ptr &gap_component,
                      EigenMatrix &A) const;

        static void fillInverseSigmaK_nm(const std::vector<IGapComponent::Ptr> &gap_components,
                                         const NeighbourList &neighbour_list,
                                         const EnergyData &energy_data,
                                         const Regularization &sigmas_inverse,
                                         EigenMatrix &A,
                                         size_t starting_row);

        static EigenMatrix choleskyDecomposition(EigenMatrix& matrix_block);
        static EigenMatrix convertToEigen(Matrix& matrix_block);
    };
}

#endif