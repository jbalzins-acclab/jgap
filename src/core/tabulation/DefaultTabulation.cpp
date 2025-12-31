#include "core/tabulation/SimpleTabulation.hpp"

#include "io/log/CurrentLogger.hpp"
#include "io/tabgap/TabGapIO.hpp"
#include "utils/AtomicNumbers.hpp"
#include "utils/BSplineTools.hpp"

namespace jgap {
    std::shared_ptr<SimpleTabulation> SimpleTabulation::fromDataNode(const DataNode &params) {
        OPTIONALLY_SET(table_filename_prefix_, params, "table_filename_prefix");

        default_params_.n2b = params.value("n2b", 20000);
        default_params_.n3bR = params.value("n3b_r", 150);
        default_params_.n3bR = params.value("n3b_angles", 150);
        default_params_.nEamDensities = params.value("n_eam_densities", 20000);

        OPTIONALLY_SET(default_params_.minDensity, params, "min_eam_density");
        OPTIONALLY_SET(default_params_.maxDensity, params, "max_eam_density");
        OPTIONALLY_SET(default_params_.cutoff2b, params, "r_max_2b");
        OPTIONALLY_SET(default_params_.cutoff3b, params, "r_max_3b");
        OPTIONALLY_SET(default_params_.rMin3b, params, "r_min_3b");
    }

    SimpleTabulation::SimpleTabulation(TabulationParams defaultParams, std::optional<std::string> tableFilenamePrefix) {

    }

    std::shared_ptr<TabGapPotential> SimpleTabulation::tabulate(const std::shared_ptr<Potential> &potential) {

        JGAP_LOG_INFO("Starting tabulation");
        auto [valueTables, splineTables] = makeSplineTables(potential, prepareParams(potential));

        JGAP_LOG_INFO("Saving tabGAP files");
        auto tableFiles = TabGapIO::write(valueTables, splineTables, table_filename_prefix_);

        return std::make_shared<TabGapPotential>(splineTables, tableFiles);
    }

    std::pair<TabulationData, TabulationData> SimpleTabulation::makeSplineTables(
        const std::shared_ptr<Potential> &potential, const TabulationParams& params) {

        JGAP_LOG_DEBUG("Starting value tabulation");
        TabulationData valueTables(params);
        potential->tabulate(valueTables);

        JGAP_LOG_DEBUG("Converting to spline coefficients");

        TabulationData splineTables;
        splineTables.isolatedEnergies = valueTables.isolatedEnergies;

        for (const auto& [pair, energies]: valueTables.pairGrids) {
            splineTables.pairGrids[pair] = BSplineTools::toSplineCoefficients(energies);
        }
        for (const auto& [triplet, energies]: valueTables.tripletGrids) {
            splineTables.tripletGrids[triplet] = BSplineTools::toSplineCoefficients(energies);
        }
        for (const auto& valuePart: valueTables.eamTabulationData) {
            auto splineEamTable = splineTables.newEamGrids();
            for (const auto& [species, embeddingEnergies]: valuePart.densityGrids) {
                splineEamTable.densityGrids[species] = BSplineTools::toSplineCoefficients(embeddingEnergies);
            }
            for (const auto& [pair, func]: valuePart.eamPairFunctionGrids) {
                splineEamTable.eamPairFunctionGrids[pair] = BSplineTools::toSplineCoefficients(func);
            }
        }

        return {valueTables, splineTables};
    }

    TabulationParams SimpleTabulation::prepareParams(const std::shared_ptr<Potential> &potential) const {
        TabulationParams params = default_params_;

        CutoffRanges cutoffs = potential->getCutoff();
        if (!params.cutoff2b.has_value()) {
            params.cutoff2b = cutoffs.twoBody;
        } else if (cutoffs.twoBody.has_value() && abs(params.cutoff2b.value() - cutoffs.twoBody.value()) > 1e-5) {
            JGAP_LOG_WARN(
                "Specified 2b cutoff={} doesn't match the potential's 2b cutoff {}",
                params.cutoff2b.value(), cutoffs.twoBody.value()
                );
        }
        if (!params.cutoff3b.has_value()) {
            params.cutoff3b = cutoffs.threeBody;
        } else if (cutoffs.threeBody.has_value()
                    && abs(params.cutoff3b.value() - cutoffs.threeBody.value()) > 1e-5) {
            JGAP_LOG_WARN(
                "Specified 3b cutoff={} doesn't match the potential's 3b cutoff={}",
                params.cutoff3b.value(), cutoffs.threeBody.value()
            );
        }
        if (!params.minDensity.has_value()) {
            params.minDensity = cutoffs.minEam;
        } else if (cutoffs.minEam.has_value() && abs(params.minDensity.value() - cutoffs.minEam.value()) > 1e-5) {
            JGAP_LOG_WARN(
                "Specified minDensity={} doesn't match the potential's minDensity={}",
                params.minDensity.value(), cutoffs.minEam.value()
            );
        }
        if (!params.maxDensity.has_value()) {
            params.maxDensity = cutoffs.maxEam;
        } else if (cutoffs.maxEam.has_value() && abs(params.maxDensity.value() - cutoffs.maxEam.value()) > 1e-5) {
            JGAP_LOG_WARN(
                "Specified maxDensity={} doesn't match the potential's maxDensity={}",
                params.maxDensity.value(), cutoffs.maxEam.value()
            );
        }
        return params;
    }
}
