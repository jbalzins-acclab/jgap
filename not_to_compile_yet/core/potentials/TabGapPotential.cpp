#include "core/potentials/TabGapPotential.hpp"

#include <utils/BSplineTools.hpp>

#include "core/cutoff/FakeCutoff.hpp"
#include "core/descriptors/2b/TabKernel2b.hpp"
#include "core/descriptors/3b/TabKernel3b.hpp"
#include "core/descriptors/eam/TabKernelEam.hpp"
#include "core/descriptors/eam/pair_functions/FailOnUsePairFunction.hpp"
#include "core/descriptors/eam/pair_functions/TabPairFunction.hpp"
#include "../../../include/data/atomic/AtomicStructure.hpp"
#include "data/PredictionData.hpp"
#include "io/tabgap/TabGapIO.hpp"

namespace jgap {
    TabGapPotential::TabGapPotential(const DataNode &params) : params_(params) {

        if (!params.contains("files")) {
            JGAP_LOG_AND_THROW("tabGAP params do not indicate files: {}", params.dump());
        }

        TabulationData splineCoefficients = TabGapIO::read(params["files"]);
        init(splineCoefficients);
    }

    TabGapPotential::TabGapPotential(TabulationData spline_coefficients, const std::vector<std::string>& files) {
        params_["files"] = files;
        init(spline_coefficients);
    }

    Predictions TabGapPotential::predict(const AtomicStructure &structure) {

        double energy = 0.0;
        for (const auto& atom: structure) {
            energy += GET_OR_DEFAULT(isolated_energies_, atom.species(), 0.0);
        }
        Predictions result = {energy, {}, {}};

        if (two_body_ != nullptr) result = result + two_body_->predict(structure);
        if (three_body_ != nullptr) result = result + three_body_->predict(structure);
        for (const auto& eam: eams_) {
            result = result + eam->predict(structure);
        }

        return result;
    }

    DataNode TabGapPotential::serialize() {
        return params_;
    }

    CutoffRanges TabGapPotential::getCutoff() {
        CutoffRanges result = {};

        if (two_body_ != nullptr) result = result + two_body_->getCutoff();
        if (three_body_ != nullptr) result = result + three_body_->getCutoff();
        for (const auto& eam: eams_) {
            result = result + eam->getCutoff();
        }

        return result;
    }

    void TabGapPotential::tabulate(TabulationData &table) {
        if (two_body_ != nullptr) two_body_->tabulate(table);
        if (three_body_ != nullptr) three_body_->tabulate(table);
        for (const auto& eam: eams_) {
            eam->tabulate(table);
        }
    }

    void TabGapPotential::init(TabulationData& spline_coefficients) {

        double cutoff2b = 0.0;
        std::vector<std::shared_ptr<TwoBodyKernel>> kernels2b;
        for (auto &[speciesPair, grid]: spline_coefficients.pairGrids) {
            cutoff2b = std::max(cutoff2b, grid.cutoff());
            kernels2b.push_back(std::make_shared<TabKernel2b>(speciesPair, std::make_shared<Grid1d>(std::move(grid))));
        }
        two_body_ = std::make_shared<TwoBodyDescriptor>(std::make_shared<FakeCutoff>(cutoff2b), kernels2b);

        double cutoff3b = 0.0;
        std::vector<std::shared_ptr<ThreeBodyKernel>> kernels3b;
        for (auto &[speciesTriplet, grid]: spline_coefficients.tripletGrids) {
            cutoff3b = std::max(cutoff3b, grid.cutoff());
            kernels3b.push_back(std::make_shared<TabKernel3b>(speciesTriplet, std::make_shared<Grid3d>(std::move(grid))));
        }
        three_body_ = std::make_shared<ThreeBodyDescriptorFinder>(std::make_shared<FakeCutoff>(cutoff3b), kernels3b);

        for (auto& eamPart: spline_coefficients.eamTabulationData) {

            std::vector<std::shared_ptr<EamKernel>> kernelsEam;
            for (auto& [species, grid]: eamPart.densityGrids) {
                kernelsEam.push_back(std::make_shared<TabKernelEam>(species, std::make_shared<Grid1d>(std::move(grid))));
            }

            std::map<ContributorReceiverSpecies, std::shared_ptr<EamPairFunction>> pairFunctionsEam;
            for (auto& [species, grid]: eamPart.eamPairFunctionGrids) {
                pairFunctionsEam[species] = std::make_shared<TabPairFunction>(std::make_shared<Grid1d>(std::move(grid)));
            }

            std::shared_ptr<EamPairFunction> failOnUse = std::make_shared<FailOnUsePairFunction>();
            eams_.push_back(std::make_shared<EamDescriptor>(kernelsEam, failOnUse, pairFunctionsEam));
        }
    }
}
