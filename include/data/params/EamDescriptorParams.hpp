#ifndef EAMDESCRIPTORPARAMS_HPP
#define EAMDESCRIPTORPARAMS_HPP

#include <map>
#include <optional>

#include "DescriptorParams.hpp"
#include "data/BasicDataTypes.hpp"

using namespace std;

namespace jgap {

    struct EamDensityCalculationParams {

        enum class PairFunction {
            FSGEN,
            POLYCUTOFF,
            COSCUTOFF
        } pairFunction;

        enum class SpeciesUse {
            BLIND,
            EAM,
            FSSYM,
            FSGEN
        } speciesUse;

        double cutoff;

        optional<double> rmin;
        optional<double> degree;
    };

    struct EamDescriptorParams : DescriptorParams {

        optional<string> name;

        enum class KernelType {
            GAUSS
        } kernelType = KernelType::GAUSS;

        enum class SparsificationMethod {
            FULL_GRID_UNIFORM, // (..)(TODO) to cutoff
            SAMPLE_SPACE_UNIFORM,
            EQUI_DENSE // ~ equal number of training points between the sparse points
        } sparsificationMethod = SparsificationMethod::FULL_GRID_UNIFORM;

        size_t nSparsePoints = 20;
        pair<optional<double>, optional<double>> sparseRange;

        // e.g. Gaussian width
        double lengthScale = 1.0;
        // prefactor before the kernel
        double energyScale = 1.0;

        optional<EamDensityCalculationParams> defaultDensityCalculationParams;
        map<pair<Species, Species>, EamDensityCalculationParams> densityCalculationParamsPerSpecies;
    };
}

#endif //EAMDESCRIPTORPARAMS_HPP
