#ifndef JGAP_QRGAPFIT_HPP
#define JGAP_QRGAPFIT_HPP

#include "data/Vector3.hpp"
#include "io/log/CurrentLogger.hpp"
#include "core/potentials/Potential.hpp"
#include "core/descriptors/Descriptor.hpp"
#include "core/potentials/GapPotential.hpp"

#include <Eigen/Dense>

#include <memory>

#include "Fit.hpp"
#include "core/matrices/regularization/RegularizationRules.hpp"

namespace jgap {

    class QRGapFit : public Fit {
    public:
        SETUP_PARSER(Fit, QRGapFit, qr_gap)

        QRGapFit(
            const std::map<std::string, std::shared_ptr<Descriptor>> &descriptors,
            const std::shared_ptr<RegularizationRules> &regularization_rules,
            double jitter = 1e-8
            );

        std::shared_ptr<Potential> fit(const std::vector<AtomicStructure>& training_data) override;

    protected:
        std::vector<double> leastSquares(Eigen::MatrixXd &A, Eigen::VectorXd &b);

    private:
        std::map<std::string, std::shared_ptr<Descriptor>> descriptors_;
        std::shared_ptr<RegularizationRules> regularization_rules_;
        double _jitter;

        // TODO: 4-element matrix has only ~12% non-zero values i.e. ~12GB
        // => try to configure sparse matrix storage + linalg
        // * either specify manually to use it, or auto-select if #elements > 2 (?)

        [[nodiscard]]
        Eigen::MatrixXd makeA(const std::vector<std::shared_ptr<Descriptor>>& descriptors,
                              const std::vector<AtomicStructure>& atomic_structures) const;

        static Eigen::VectorXd makeB(const std::vector<std::shared_ptr<Descriptor>>& descriptors,
                                     const std::vector<AtomicStructure>& atomic_structures);

        void fillU_mm(size_t startingRow, size_t startingCol, Descriptor &descriptor, Eigen::MatrixXd &A) const;

        static void fillInverseSigmaK_nm(
            const std::vector<std::shared_ptr<Descriptor>> &descriptors,
            const AtomicStructure &atomic_structure,
            Eigen::MatrixXd &A,
            size_t starting_row
        );

        static Eigen::MatrixXd choleskyDecomposition(Eigen::MatrixXd& matrix_block);
        static Eigen::MatrixXd convertToEigen(MatrixBlock& matrix_block);
    };
}

#endif