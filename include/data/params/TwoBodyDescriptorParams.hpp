#ifndef TWOBODYDESCRIPTORPARAMS_HPP
#define TWOBODYDESCRIPTORPARAMS_HPP

#include "DescriptorParams.hpp"
#include "data/BasicDataTypes.hpp"

#include <map>
#include <string>

using namespace std;

namespace jgap {
    struct TwoBodyDescriptorParams : DescriptorParams {

        optional<string> name;

        enum class KernelType { GAUSS } kernelType = KernelType::GAUSS;

        enum class SparsificationMethod {
            FULL_GRID_UNIFORM,
            SAMPLE_SPACE_UNIFORM,
            EQUIDENSE // ~ equal number of training points between the sparse points
            } sparsificationMethod = SparsificationMethod::FULL_GRID_UNIFORM;

        double cutoff = 5; // descriptor-specific cutoff
        double cutoffTransitionWidth = 0.5;

        // TODO(not a necessity, but might be handy e.g. for a case when species aren't atoms(?)):
        //     implement ComplexTwoBodyDescriptorParams with lengthScale/energyScale per species pair.
        // prefactor before the kernel; TODO: is it not absorbed into coefficients?
        double energyScale = 10;
        // e.g. Gaussian's width
        double lengthScale = 1;

        size_t nSparsePointsPerSpeciesPair = 20;

        pair<optional<double>, optional<double>> sparseRange;

        optional<vector<SpeciesPair>> speciesPairs;
    };
}

#endif
