#ifndef CORE_HPP
#define CORE_HPP

#include "data/BasicDataTypes.hpp"
#include "io/log/Logger.hpp"
#include "potentials/Potential.hpp"
#include "core/descriptors/Descriptor.hpp"
#include "potentials/JgapPotential.hpp"
#include "matrices/sigmas/SigmaRules.hpp"

#include <Eigen/Dense>

#include <memory>
#include <concepts>


using namespace std;

namespace jgap {

    class InRamJgapFit {
    public:
        InRamJgapFit(const vector<shared_ptr<Descriptor>> &descriptors,
                     const vector<shared_ptr<Potential>> &externalPotentials = {},
                     const shared_ptr<SigmaRules>& sigmaRules = nullptr,
                     double jitter = 1e-8);

        ~InRamJgapFit() = default;

        [[nodiscard]]
        shared_ptr<JgapPotential> fitGAP(const vector<AtomicStructure>& trainingData) const;

    private:
        vector<shared_ptr<Potential>> _externalPotentials;
        vector<shared_ptr<Descriptor>> _descriptors;
        shared_ptr<SigmaRules> _sigmaRules;
        double _jitter;

        [[nodiscard]]
        vector<AtomicStructure> subtractExternalContributions(const vector<AtomicStructure> &originalData) const;

        Eigen::MatrixXd makeA(const vector<shared_ptr<Descriptor>>& descriptors,
                              vector<AtomicStructure>& atomicStructures) const;

        static Eigen::VectorXd makeB(const vector<shared_ptr<Descriptor>>& descriptors,
                              const vector<AtomicStructure>& atomicStructures) ;

        void fillU_mm(size_t startingRow, size_t startingCol, Descriptor &descriptor, Eigen::MatrixXd &A) const;

        static void fillInverseRootSigmaK_nm(
            const vector<shared_ptr<Descriptor>> &descriptors,
            const AtomicStructure &atomicStructure,
            Eigen::MatrixXd &A,
            size_t startingRow
        );

        static Eigen::MatrixXd choleskyDecomposition(Eigen::MatrixXd& matrixBlock);
        static Eigen::MatrixXd convertToEigen(MatrixBlock& matrixBlock);
    };
}

#endif