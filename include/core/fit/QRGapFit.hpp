#ifndef QRGAPFIT_HPP
#define QRGAPFIT_HPP

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
        static constexpr std::string TYPE = "qr_gap";
        static std::shared_ptr<QRGapFit> fromJson(const nlohmann::json& params);

        QRGapFit(
            const std::map<std::string, std::shared_ptr<Descriptor>> &descriptors,
            const std::shared_ptr<RegularizationRules> &regularizationRules,
            double jitter = 1e-8
            );

        std::shared_ptr<Potential> fit(const std::vector<AtomicStructure>& trainingData) override;

    protected:
        std::vector<double> leastSquares(Eigen::MatrixXd &A, Eigen::VectorXd &b);

    private:
        std::map<std::string, std::shared_ptr<Descriptor>> _descriptors;
        std::shared_ptr<RegularizationRules> _regularizationRules;
        double _jitter;

        // TODO: 4-element matrix has only ~12% non-zero values i.e. ~12GB
        // => try to configure sparse matrix storage + linalg
        // * either specify manually to use it, or auto-select if #elements > 2 (?)

        [[nodiscard]]
        Eigen::MatrixXd makeA(const std::vector<std::shared_ptr<Descriptor>>& descriptors,
                              const std::vector<AtomicStructure>& atomicStructures) const;

        static Eigen::VectorXd makeB(const std::vector<std::shared_ptr<Descriptor>>& descriptors,
                                     const std::vector<AtomicStructure>& atomicStructures);

        void fillU_mm(size_t startingRow, size_t startingCol, Descriptor &descriptor, Eigen::MatrixXd &A) const;

        static void fillInverseSigmaK_nm(
            const std::vector<std::shared_ptr<Descriptor>> &descriptors,
            const AtomicStructure &atomicStructure,
            Eigen::MatrixXd &A,
            size_t startingRow
        );

        static Eigen::MatrixXd choleskyDecomposition(Eigen::MatrixXd& matrixBlock);
        static Eigen::MatrixXd convertToEigen(MatrixBlock& matrixBlock);
    };

    REGISTER_PARSER(Fit, QRGapFit)
}

#endif