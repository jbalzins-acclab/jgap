#ifndef EAMDESCRIPTOR_HPP
#define EAMDESCRIPTOR_HPP

#include <nlohmann/json.hpp>

#include "core/descriptors/Descriptor.hpp"
#include "core/kernels/Kernel.hpp"
#include "data/kernels/EamKernelIndex.hpp"
#include "data/params/EamDescriptorParams.hpp"
#include "eam/EamDensityCalculator.hpp"
#include "memory/MatrixBlock.hpp"

namespace jgap {

    class EamDescriptor : public Descriptor {
    public:
        explicit EamDescriptor(const EamDescriptorParams& params);

        ~EamDescriptor() override = default;

        void setSparsePoints(const vector<AtomicStructure> &fromData) override;

        void setSparsePoints(vector<double> sparsePoints) {
            _sparsePoints = std::move(sparsePoints);
        }

        double getCutoff() override { return _densityCalculator->getCutoff(); }

        size_t nSparsePoints() override;

        vector<Covariance> covariate(const AtomicStructure &atomicStructure) override;
        vector<pair<size_t, shared_ptr<MatrixBlock>>> selfCovariate() override;

        nlohmann::json serialize() override;
    private:
        string _name;
        EamDescriptorParams _params;
        shared_ptr<Kernel<double, EamKernelIndex>> _kernel;
        shared_ptr<EamDensityCalculator> _densityCalculator;

        vector<double> _sparsePoints;

        void sparsifyTrueUniform(const vector<AtomicStructure> &fromData);
        void sparsifyQuipUniform(const vector<AtomicStructure> &fromData);
    };
}



#endif //EAMDESCRIPTOR_HPP
