#ifndef INRAMJGAPFIT_HPP
#define INRAMJGAPFIT_HPP

#include "data/BasicDataTypes.hpp"
#include "io/log/CurrentLogger.hpp"
#include "core/potentials/Potential.hpp"
#include "core/descriptors/Descriptor.hpp"
#include "core/potentials/JgapPotential.hpp"
#include "core/matrices/sigmas/SigmaRules.hpp"

#include <Eigen/Dense>

#include <memory>
#include <concepts>

#include "Fit.hpp"

using namespace std;

namespace jgap {

    class InRamJgapFit : public Fit {
    public:
        ~InRamJgapFit() override = default;

        explicit InRamJgapFit(const nlohmann::json& params);
        string getType() override { return "in_ram_jgap"; }

        [[nodiscard]]
        shared_ptr<Potential> fit(const vector<AtomicStructure>& trainingData) override;

    protected:
        vector<double> leastSquares(Eigen::MatrixXd &A, Eigen::VectorXd &b);

    private:
        map<string, shared_ptr<Descriptor>> _descriptors;
        shared_ptr<SigmaRules> _sigmaRules;
        double _jitter;

        [[nodiscard]]
        Eigen::MatrixXd makeA(const vector<shared_ptr<Descriptor>>& descriptors,
                              const vector<AtomicStructure>& atomicStructures) const;

        static Eigen::VectorXd makeB(const vector<shared_ptr<Descriptor>>& descriptors,
                                     const vector<AtomicStructure>& atomicStructures);

        void fillU_mm(size_t startingRow, size_t startingCol, Descriptor &descriptor, Eigen::MatrixXd &A) const;

        static void fillInverseSigmaK_nm(
            const vector<shared_ptr<Descriptor>> &descriptors,
            const AtomicStructure &atomicStructure,
            Eigen::MatrixXd &A,
            size_t startingRow
        );

        static Eigen::MatrixXd choleskyDecomposition(Eigen::MatrixXd& matrixBlock);
        static Eigen::MatrixXd convertToEigen(MatrixBlock& matrixBlock);
    };

    REGISTER_PARSER("in_ram_jgap", Fit, InRamJgapFit)
}

#endif