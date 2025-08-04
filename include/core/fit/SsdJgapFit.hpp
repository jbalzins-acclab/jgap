#ifndef SSDJGAPFIT_HPP
#define SSDJGAPFIT_HPP
#ifndef _WIN32

#include "data/BasicDataTypes.hpp"
#include "io/log/CurrentLogger.hpp"
#include "core/potentials/Potential.hpp"
#include "core/descriptors/Descriptor.hpp"
#include "core/potentials/JgapPotential.hpp"
#include "core/matrices/sigmas/RegularizationRules.hpp"

#include <Eigen/Dense>

#include <memory>
#include <concepts>

#include "core/fit/Fit.hpp"

using namespace std;

namespace jgap {

    class SsdJgapFit : public Fit {
    public:
        ~SsdJgapFit() override = default;

        explicit SsdJgapFit(const nlohmann::json& params);
        string getType() override { return "ssd_jgap"; }

        shared_ptr<Potential> fit(const vector<AtomicStructure>& trainingData) override;

    private:
        string _ssdTmpDir;

        map<string, shared_ptr<Descriptor>> _descriptors;
        shared_ptr<RegularizationRules> _regularizationRules;
        double _jitter;

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

    REGISTER_PARSER("ssd_jgap", Fit, SsdJgapFit)
}

#endif
#endif