#include <core/InRamJgapFit.hpp>
#include <utils/Utils.hpp>
#include <iostream>
#include <core/neighbours/NeighbourFinder.hpp>

#include <algorithm>
#include <fstream>
#include <nlohmann/json.hpp>

#include "core/descriptors/EamDescriptor.hpp"
#include "core/descriptors/ThreeBodyDescriptor.hpp"
#include "core/descriptors/TwoBodyDescriptor.hpp"
#include "core/potentials/IsolatedAtomPotential.hpp"
#include "core/potentials/ZblPotential.hpp"
#include "data/params/EamDescriptorParams.hpp"
#include "data/params/ThreeBodyDescriptorParams.hpp"
#include "io/log/StdoutLogger.hpp"
#include <format>
#include <Eigen/Core>

#include "core/matrices/sigmas/SimpleSigmaRules.hpp"
using namespace std;

int main() {

    Eigen::setNbThreads(6);

    const auto params2b = jgap::TwoBodyDescriptorParams{
        .cutoff = 4.5,
        .cutoffTransitionWidth = 1.0,
        .kernelType = jgap::TwoBodyDescriptorParams::KernelType::GAUSS,
        .sparsificationMethod = jgap::TwoBodyDescriptorParams::SparsificationMethod::SAMPLE_SPACE_UNIFORM,
        .nSparsePointsPerSpeciesPair = 20,
        .energyScale = 10.0,
        .lengthScale = 1.0
    };
    shared_ptr<jgap::Descriptor> desc2b = make_shared<jgap::TwoBodyDescriptor>(params2b);

    constexpr auto densityParams = jgap::EamDensityCalculationParams{
        .pairFunction = jgap::EamDensityCalculationParams::PairFunction::POLYCUTOFF,
        .cutoff = 4.5,
        .rmin = 0,
        .speciesUse = jgap::EamDensityCalculationParams::SpeciesUse::BLIND
    };

    const auto paramsEAM = jgap::EamDescriptorParams{
        .kernelType = jgap::EamDescriptorParams::KernelType::GAUSS,
        .sparsificationMethod = jgap::EamDescriptorParams::SparsificationMethod::SAMPLE_SPACE_UNIFORM,
        .nSparsePoints = 20,
        .energyScale = 1.0,
        .lengthScale = 1.0,
        .defaultDensityCalculationParams = densityParams
    };
    shared_ptr<jgap::Descriptor> descEAM = make_shared<jgap::EamDescriptor>(paramsEAM);

    const auto params3b = jgap::ThreeBodyDescriptorParams{
        .kernelType = jgap::ThreeBodyDescriptorParams::KernelType::GAUSS,
        .sparsificationMethod = jgap::ThreeBodyDescriptorParams::SparsificationMethod::SAMPLE_SPACE_UNIFORM,
        .nSparsePointsPerSpecies = 500,
        .energyScale = 1.0,
        .lengthScale = 1.0,
        .cutoff = 3.7,
        .cutoffTransitionWidth = 0.6,
    };
    shared_ptr<jgap::Descriptor> desc3b = make_shared<jgap::ThreeBodyDescriptor>(params3b);

    string fn = "resources/xyz-samples/FeOnly.xyz";
    auto structsAll = jgap::readXyz(fn);
    jgap::NeighbourFinder::findNeighbours(structsAll, 5.0);

    vector<jgap::AtomicStructure> structs;
    for (int i = 0; i < structsAll.size(); ++i) {
        if (structsAll[i].configType.value_or("-") == "isolated_atom") {
            structs.push_back(structsAll[i]);
            // structs.back().energySigma = 0.00001;//??
        } else if (true) {//(67 < i && i < 437) {
            structs.push_back(structsAll[i]);
        }
    }

    shared_ptr<jgap::Potential> isolatedAtomPotential = make_shared<jgap::IsolatedAtomPotential>(structsAll);

    ifstream file("resources/dmol-screening-fit/dmol-fit.json");
    nlohmann::json dmolData{};
    file >> dmolData;
    shared_ptr<jgap::Potential> zblPotential = make_shared<jgap::ZblPotential>(dmolData);

    jgap::InRamJgapFit fit = jgap::InRamJgapFit(
        {desc2b, descEAM, desc3b},
        {isolatedAtomPotential, zblPotential},
        make_shared<jgap::SimpleSigmaRules>(0.001, 0.04, 10, 10));
    //jgap::InRamJgapFit fit = jgap::InRamJgapFit(descs, {isolatedAtomPotential, zblPotential});
    auto res = fit.fitGAP(structs);

    ofstream fOut(format("iron{}.json", "-2b-eam-3b"));
    fOut << res->serialize().dump() << endl;
    fOut.close();

    int ii = 15;
    auto pred = res->predict(structsAll[ii])
                               + zblPotential -> predict(structsAll[ii])
                                + isolatedAtomPotential->predict(structsAll[ii]);
    jgap::Logger::logger->info(format(
        "Prediction test: was {} now {}",
        structsAll[ii].energy.value(),
        pred.energy.value()
        ));
    for (size_t i = 0; i < structsAll[ii].atoms.size(); ++i) {
        jgap::Logger::logger->info(format(
            "f[{}]={}quip|{}jgap", i,
            structsAll[ii].atoms[i].force.value().toString(),
            pred.forces.value()[i].toString()
            ));
    }

    //auto fit = jgap::JgapFit(, nullptr, nullptr, {});
    // jgap::Logger::logger->info(res->serialize().dump());

    return 0;
}
