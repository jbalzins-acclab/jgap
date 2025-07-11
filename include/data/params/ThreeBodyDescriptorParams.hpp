#ifndef THREEBODYDESCRIPTORPARAMS_HPP
#define THREEBODYDESCRIPTORPARAMS_HPP

#include <map>

#include "DescriptorParams.hpp"
#include "data/BasicDataTypes.hpp"

using namespace std;

namespace jgap {
    class ThreeBodyDescriptorParams : public DescriptorParams {
    public:

        optional<string> name;

        enum class KernelType { GAUSS } kernelType = KernelType::GAUSS;

        enum class SparsificationMethod {
            FULL_GRID_UNIFORM,
            SAMPLE_SPACE_UNIFORM
            // TODO: Clustering?
            } sparsificationMethod = SparsificationMethod::FULL_GRID_UNIFORM;

        double cutoff = 5;
        double cutoffTransitionWidth = 0.5;

        // e.g. Gaussian's width:
        double lengthScale = 1;
        // prefactor before the kernel; TODO: is it not absorbed into coefficients?
        double energyScale = 2;

        // TODO: sub-class?
        array<array<double, 2/*from, to*/>, 3/*spec in each direction*/> sparseRanges;
        optional<array<size_t, 3>> nSparsePointsPerSpeciesPerDirection;
        optional<size_t> nSparsePointsPerSpecies;

        optional<vector<SpeciesTriplet>> speciesTriplets;
    };
}

#endif
