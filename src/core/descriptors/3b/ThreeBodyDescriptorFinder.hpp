#ifndef THREEBODYDESCRIPTOR_HPP
#define THREEBODYDESCRIPTOR_HPP

#include <queue>
#include <utility>

#include "core/descriptors/Descriptor.hpp"
#include "../../kernels/3b/ThreeBodySE.hpp"
#include "core/descriptors/DescriptorFinder.hpp"
#include "core/kernels/KernelCollection.hpp"
#include "core/sparsification/Sparsifier.hpp"
#include "core/tabulation/Tabulatable.hpp"
#include "data/descriptors/kernels/ThreeBodyIndex.hpp"

#include "io/parse/ParserRegistry.hpp"
#include "memory/MatrixBlock.hpp"

namespace jgap {

    class ThreeBodyDescriptorFinder/*, Serializable, Tabulatable*/ {
    public:
        //SETUP_PARSER_AND_SERIALIZATION(ThreeBodyDescriptorFinder, ThreeBodyDescriptorFinder, 3b)

        using DescriptorsFiltered = std::map<EncodedSpeciesSets, >

        DescriptorsFiltered find(const AtomicStructure &atomic_structure, double cutoff) override;

        Vector3 toInvariantTriplet(double r01, double r02, double r12) const;
        std::array<Vector3, 3> invariantTripletGradients(double r01, double r02) const;

        //void tabulate(TabulationData &table) override;
    };

    class ThreeBodyKernelCollection : public KernelCollection, Tabulatable, Serializable {
    public:
        std::vector<Predictions> covariate(const AtomicStructure &atomic_structure) override;

        std::vector<std::shared_ptr<MatrixBlock>> selfCovariate() override;

        CutoffRanges getCutoff() override { return {.threeBody = cutoff_function_->getCutoff()}; }

        std::vector<std::shared_ptr<IKernel>> getKernels() override;

    private:
        std::shared_ptr<CutoffFunction> cutoff_function_;
        std::shared_ptr<ThreeBodyDescriptorFinder> descriptor_finder_;
        std::map<EncodedSpeciesSets, std::vector<std::shared_ptr<ThreeBodyKernel>>> kernels_per_species_triplet_;
    };
}

#endif
