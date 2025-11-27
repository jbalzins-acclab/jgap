#include "core/potentials/TabGapPotential.hpp"

#include <utils/BSplineTools.hpp>

#include "core/cutoff/FakeCutoff.hpp"
#include "core/descriptors/2b/TabKernel2b.hpp"
#include "core/descriptors/3b/TabKernel3b.hpp"
#include "core/descriptors/eam/TabKernelEam.hpp"
#include "core/descriptors/eam/pair_functions/FailOnUsePairFunction.hpp"
#include "core/descriptors/eam/pair_functions/TabPairFunction.hpp"
#include "data/AtomicStructure.hpp"
#include "data/PredictionData.hpp"
#include "io/tabgap/TabGapIO.hpp"

namespace jgap {
    TabGapPotential::TabGapPotential(const nlohmann::json &params) : _params(params) {

        if (!params.contains("files")) {
            JGAP_LOG_AND_THROW("tabGAP params do not indicate files: {}", params.dump());
        }

        TabulationData splineCoefficients = TabGapIO::read(params["files"]);
        init(splineCoefficients);
    }

    TabGapPotential::TabGapPotential(TabulationData splineCoefficients, const vector<string>& files) {
        _params["files"] = files;
        init(splineCoefficients);
    }

    Predictions TabGapPotential::predict(const AtomicStructure &structure) {

        double energy = 0.0;
        for (const auto& atom: structure) {
            energy += GET_OR_DEFAULT(_isolatedEnergies, atom.species(), 0.0);
        }
        Predictions result = {energy, {}, {}};

        if (_2b != nullptr) result = result + _2b->predict(structure);
        if (_3b != nullptr) result = result + _3b->predict(structure);
        for (const auto& eam: _eams) {
            result = result + eam->predict(structure);
        }

        return result;
    }

    nlohmann::json TabGapPotential::serialize() {
        return _params;
    }

    CutoffRanges TabGapPotential::getCutoff() {
        CutoffRanges result = {};

        if (_2b != nullptr) result = result + _2b->getCutoff();
        if (_3b != nullptr) result = result + _3b->getCutoff();
        for (const auto& eam: _eams) {
            result = result + eam->getCutoff();
        }

        return result;
    }

    void TabGapPotential::tabulate(TabulationData &table) {
        if (_2b != nullptr) _2b->tabulate(table);
        if (_3b != nullptr) _3b->tabulate(table);
        for (const auto& eam: _eams) {
            eam->tabulate(table);
        }
    }

    void TabGapPotential::init(TabulationData& splineCoefficients) {

        double cutoff2b = 0.0;
        vector<shared_ptr<TwoBodyKernel>> kernels2b;
        for (auto &[speciesPair, grid]: splineCoefficients.pairGrids) {
            cutoff2b = max(cutoff2b, grid.cutoff());
            kernels2b.push_back(make_shared<TabKernel2b>(speciesPair, make_shared<Grid1d>(std::move(grid))));
        }
        _2b = make_shared<TwoBodyDescriptor>(make_shared<FakeCutoff>(cutoff2b), kernels2b);

        double cutoff3b = 0.0;
        vector<shared_ptr<ThreeBodyKernel>> kernels3b;
        for (auto &[speciesTriplet, grid]: splineCoefficients.tripletGrids) {
            cutoff3b = max(cutoff3b, grid.cutoff());
            kernels3b.push_back(make_shared<TabKernel3b>(speciesTriplet, make_shared<Grid3d>(std::move(grid))));
        }
        _3b = make_shared<ThreeBodyDescriptor>(make_shared<FakeCutoff>(cutoff3b), kernels3b);

        for (auto& eamPart: splineCoefficients.eamTabulationData) {

            vector<shared_ptr<EamKernel>> kernelsEam;
            for (auto& [species, grid]: eamPart.densityGrids) {
                kernelsEam.push_back(make_shared<TabKernelEam>(species, make_shared<Grid1d>(std::move(grid))));
            }

            map<ContributorReceiverSpecies, shared_ptr<EamPairFunction>> pairFunctionsEam;
            for (auto& [species, grid]: eamPart.eamPairFunctionGrids) {
                pairFunctionsEam[species] = make_shared<TabPairFunction>(make_shared<Grid1d>(std::move(grid)));
            }

            shared_ptr<EamPairFunction> failOnUse = make_shared<FailOnUsePairFunction>();
            _eams.push_back(make_shared<EamDescriptor>(kernelsEam, failOnUse, pairFunctionsEam));
        }
    }
}
