#include "io/convert/QuipXmlConverter.hpp"

#include "core/cutoff/PerriotPolynomialCutoff.hpp"
#include "core/potentials/CompositePotential.hpp"
#include "core/potentials/IsolatedAtomPotential.hpp"
#include "core/potentials/SplinePairPotential.hpp"
#include "utils/AtomicNumbers.hpp"

#include <fstream>

#include "core/cutoff/CosCutoff.hpp"
#include "core/descriptors/eam/pair_functions/CoscutoffPairFunction.hpp"
#include "core/descriptors/eam/pair_functions/FSGenPairFunction.hpp"
#include "core/descriptors/eam/pair_functions/PolycutoffPairFunction.hpp"

namespace jgap {
    std::shared_ptr<Potential> QuipXmlConverter::transform(const pugi::xml_node quipPotential) {

        if (quipPotential.child("pairpot").empty()) {
            if (quipPotential.child("GAP_params").empty()) {
                JGAP_LOG_ERROR(
                    "Neither 'pairpot' nor 'GAP_params' were found in the xml node",
                    true
                    );
            }
            return transformGapParams(quipPotential.child("GAP_params"));
        }
        if (quipPotential.child("GAP_params").empty()) {
            return transformPairpot(quipPotential.child("pairpot"));
        }

        auto pairpot = transformPairpot(quipPotential.child("pairpot"));
        auto quipGap = transformGapParams(quipPotential.child("GAP_params"));
        std::map<std::string, std::shared_ptr<Potential>> potentials = {
            {"quip_pairpot", pairpot},
            {"quip_gap", quipGap},
        };

        return std::make_shared<CompositePotential>(potentials);
    }

    std::shared_ptr<Potential> QuipXmlConverter::transformPairpot(const pugi::xml_node quipPairpot) {
        if (quipPairpot.child("Potential").attribute("init_args").as_string() != std::string("IP Glue")) {
            JGAP_LOG_WARN("Strange 'init_args'");
        }
        if (!quipPairpot.child("Glue_params")) {
            JGAP_LOG_ERROR("Non 'Glue' pairpot", true);
        }

        std::map<std::string, Species> typeToSpecies = {};
        for (pugi::xml_node perTypeNode: quipPairpot.child("Glue_params").children("per_type_data")) {
            const size_t atomicNumber = perTypeNode.attribute("atomic_num").as_uint();
            typeToSpecies[perTypeNode.attribute("type").as_string()] = CHEM_SYMBOLS[atomicNumber];
        }

        std::map<SpeciesPair, std::pair<std::vector<double>, std::vector<double>>> points;
        for (pugi::xml_node perPairNode: quipPairpot.child("Glue_params").children("per_pair_data")) {
            Species species1 = typeToSpecies[perPairNode.attribute("type1").as_string()];
            Species species2 = typeToSpecies[perPairNode.attribute("type2").as_string()];

            std::vector<double> r, E;
            for (pugi::xml_node pointNode: perPairNode.child("potential_pair").children("point")) {
                r.push_back(pointNode.attribute("r").as_double());
                E.push_back(pointNode.attribute("E").as_double());
            }

            points[SpeciesPair(species1, species2)] = std::pair{r, E};
        }

        return std::make_shared<SplinePairPotential>(std::map(points));
    }

    std::shared_ptr<Potential> QuipXmlConverter::transformGapParams(const pugi::xml_node quipGapParams) {

        std::map<std::string, std::shared_ptr<Potential>> potentials;
        if (!quipGapParams.child("GAP_data").empty()) {
            potentials["isolated_atom"] = transformIsolatedAtomParams(quipGapParams.child("GAP_data"));
        }
        if (!quipGapParams.child("gpSparse").empty()) {
            potentials["GAP"] = transformSparseData(quipGapParams.child("gpSparse"));
        }

        return std::make_shared<CompositePotential>(potentials);
    }

    bool QuipXmlConverter::QuipDescriptorData::operator==(const QuipDescriptorData &other) const {
        return std::tie(type, delta, theta, cutoff, rMin, cutoffTransitionWidth, pairFunction, order, mode)
               == std::tie(other.type, other.delta, other.theta, other.cutoff,
                           other.rMin, other.cutoffTransitionWidth, other.pairFunction, other.order, other.mode);
    }

    bool QuipXmlConverter::QuipDescriptorData::operator<(const QuipDescriptorData &other) const {
        return std::tie(type, delta, theta, cutoff, rMin, cutoffTransitionWidth, pairFunction, order, mode)
               < std::tie(other.type, other.delta, other.theta, other.cutoff,
                          other.rMin, other.cutoffTransitionWidth, other.pairFunction, other.order, other.mode);
    }

    std::shared_ptr<IsolatedAtomPotential> QuipXmlConverter::transformIsolatedAtomParams(
                                                        const pugi::xml_node quipIsolatedAtomParams) {
        std::map<Species, double> isolatedAtomEnergies;
        for (pugi::xml_node isolatedAtomNode: quipIsolatedAtomParams.children("e0")) {
            if (isolatedAtomNode.attribute("value").as_double() != 0.0) {
                const size_t atomicNumber = isolatedAtomNode.attribute("Z").as_uint();
                isolatedAtomEnergies[CHEM_SYMBOLS[atomicNumber]] = isolatedAtomNode.attribute("value").as_double();
            }
        }

        return std::make_shared<IsolatedAtomPotential>(isolatedAtomEnergies, false);
    }

    std::shared_ptr<GapPotential> QuipXmlConverter::transformSparseData(pugi::xml_node quipSparseData) {

        std::map<QuipDescriptorData, std::vector<pugi::xml_node>> nodesBySimilarity;
        for (pugi::xml_node sparseNode: quipSparseData.children("gpCoordinates")) {

            std::string descriptorParamString = sparseNode.child("descriptor").first_child().value();

            std::string type;
            if (descriptorParamString.contains("eam_density")) {
                type = "eam_density";
            } else if (descriptorParamString.contains("distance_2b")) {
                type = "distance_2b";
            } else if (descriptorParamString.contains("angle_3b")) {
                type = "angle_3b";
            } else {
                JGAP_LOG_AND_THROW("Unknown descriptor type {}", descriptorParamString);
            }

            if (!descriptorParamString.contains("covariance_type=ard_se")) {
                JGAP_LOG_AND_THROW("covariance_type must be ard_se: {}", descriptorParamString);
            }

            // TODO: make it pretty
            auto deltaStartIdx = descriptorParamString.find("delta=") + std::string("delta=").size();
            auto deltaEndIdx = descriptorParamString.find(' ', deltaStartIdx);
            std::string deltaStr = descriptorParamString.substr(deltaStartIdx, deltaEndIdx - deltaStartIdx);
            double delta = stod(deltaStr);

            auto thetaStartIdx = descriptorParamString.find("theta_uniform=") + std::string("theta_uniform=").size();
            auto thetaEndIdx = descriptorParamString.find(' ', thetaStartIdx);
            std::string thetaStr = descriptorParamString.substr(thetaStartIdx, thetaEndIdx - thetaStartIdx);
            double theta = stod(thetaStr);

            auto cutoffStartIdx = descriptorParamString.find("cutoff=") + std::string("cutoff=").size();
            auto cutoffEndIdx = descriptorParamString.find(' ', cutoffStartIdx);
            std::string cutoffStr = descriptorParamString.substr(cutoffStartIdx, cutoffEndIdx - cutoffStartIdx);
            double cutoff = stod(cutoffStr);

            std::optional<double> rMin;
            if (descriptorParamString.contains("rmin=")) {
                auto rMinStartIdx = descriptorParamString.find("rmin=") + std::string("rmin=").size();
                auto rMinEndIdx = descriptorParamString.find(' ', rMinStartIdx);
                std::string rMinStr = descriptorParamString.substr(rMinStartIdx, rMinEndIdx - rMinStartIdx);
                rMin = stod(rMinStr);
            }

            std::optional<double> cutoffTransitionWidth;
            if (descriptorParamString.contains("cutoff_transition_width=")) {
                auto cutoffTWStartIdx = descriptorParamString.find("cutoff_transition_width=")
                                                        + std::string("cutoff_transition_width=").size();
                auto cutoffTWEndIdx = descriptorParamString.find(' ', cutoffTWStartIdx);
                std::string cutoffTwStr = descriptorParamString.substr(cutoffTWStartIdx, cutoffTWEndIdx - cutoffTWStartIdx);
                cutoffTransitionWidth = stod(cutoffTwStr);
            }

            std::optional<std::string> pairFunction;
            if (descriptorParamString.contains("pair_function=")) {
                auto pfStartIdx = descriptorParamString.find("pair_function=")
                                                + std::string("pair_function=").size();
                auto pfEndIdx = descriptorParamString.find(' ', pfStartIdx);
                pairFunction = descriptorParamString.substr(pfStartIdx, pfEndIdx - pfStartIdx);
            }

            std::optional<std::string> mode;
            if (descriptorParamString.contains("mode=")) {
                auto modeStartIdx = descriptorParamString.find("mode=")
                                                + std::string("mode=").size();
                auto modeEndIdx = descriptorParamString.find(' ', modeStartIdx);
                mode = descriptorParamString.substr(modeStartIdx, modeEndIdx - modeStartIdx);
            }

            std::optional<double> order;
            if (descriptorParamString.contains("order=")) {
                auto orderStartIdx = descriptorParamString.find("order=")
                                                + std::string("order=").size();
                auto orderEndIdx = descriptorParamString.find(' ', orderStartIdx);
                order = stod(descriptorParamString.substr(orderStartIdx, orderEndIdx - orderStartIdx));
            }

            QuipDescriptorData mainData = {
                .type = type,
                .delta = delta,
                .theta = theta,
                .cutoff = cutoff,
                .rMin = rMin,
                .cutoffTransitionWidth = cutoffTransitionWidth,
                .pairFunction = pairFunction,
                .order = order,
                .mode = mode,
            };

            if (!nodesBySimilarity.contains(mainData)) {
                nodesBySimilarity[mainData] = {};
            }
            nodesBySimilarity[mainData].push_back(sparseNode);
        }

        std::map<std::string, std::shared_ptr<Descriptor>> descriptors;
        size_t cnt = 0;
        for (const auto &[mainData, nodes]: nodesBySimilarity) {
            if (mainData.type == "distance_2b") {
                descriptors[std::to_string(cnt++)] = transformDistance2b(mainData, nodes);
            } else if (mainData.type == "angle_3b") {
                descriptors[std::to_string(cnt++)] = transformAngle3b(mainData, nodes);
            } else if (mainData.type == "eam_density") {
                descriptors[std::to_string(cnt++)] = transformEam(mainData, nodes);
            } else {
                JGAP_LOG_ERROR("Unknown descriptor type: " + mainData.type, true);
            }
        }

        return std::make_shared<GapPotential>(descriptors);
    }

    std::shared_ptr<TwoBodyDescriptor> QuipXmlConverter::transformDistance2b(
        QuipDescriptorData mainData, std::vector<pugi::xml_node> distance2bNodes) {

        double cutoffTransitionWidth = mainData.cutoffTransitionWidth.value_or(0.5);
        if (mainData.rMin.has_value()) cutoffTransitionWidth = mainData.cutoff - mainData.rMin.value();

        std::shared_ptr<CutoffFunction> cutoffFunction
            = std::make_shared<CosCutoff>(mainData.cutoff, cutoffTransitionWidth);
        std::vector<std::shared_ptr<TwoBodyKernel>> kernels;

        for (pugi::xml_node distanceNode: distance2bNodes) {
            std::string descriptorParamString = distanceNode.child("descriptor").first_child().value();

            auto z1StartIdx = descriptorParamString.find("Z1=") + std::string("Z1=").size();
            auto z1EndIdx = descriptorParamString.find(' ', z1StartIdx);
            Species species1 = CHEM_SYMBOLS[stoi(descriptorParamString.substr(z1StartIdx, z1EndIdx - z1StartIdx))];

            auto z2StartIdx = descriptorParamString.find("Z2=") + std::string("Z2=").size();
            auto z2EndIdx = descriptorParamString.find(' ', z2StartIdx);
            Species species2 = CHEM_SYMBOLS[stoi(descriptorParamString.substr(z2StartIdx, z2EndIdx - z2StartIdx))];

            SpeciesPair sp{species1, species2};

            // coeffs
            double r, coeff;
            std::ifstream fin(distanceNode.attribute("sparseX_filename").as_string());
            for (pugi::xml_node pt: distanceNode.children("sparseX")) {
                fin >> r;
                coeff = pt.attribute("alpha").as_double();
                kernels.push_back(std::make_shared<TwoBodySE>(sp, mainData.delta, mainData.theta, r, coeff));
            }
            fin.close();
        }

        return std::make_shared<TwoBodyDescriptor>(cutoffFunction, kernels);
    }

    std::shared_ptr<ThreeBodyDescriptor> QuipXmlConverter::transformAngle3b(QuipDescriptorData mainData,
                                                                       std::vector<pugi::xml_node> angle3bNodes) {

        double rMin = mainData.cutoff - mainData.cutoffTransitionWidth.value_or(0.5);
        if (mainData.rMin.has_value()) rMin = mainData.rMin.value();

        std::shared_ptr<CutoffFunction> cutoffFunction = std::make_shared<CosCutoff>(mainData.cutoff, rMin);
        std::vector<std::shared_ptr<ThreeBodyKernel>> kernels;

        std::map<SpeciesTriplet, std::pair<std::vector<Vector3>, std::vector<double>>> pointsAndCoefficients;
        for (pugi::xml_node distanceNode: angle3bNodes) {
            std::string descriptorParamString = distanceNode.child("descriptor").first_child().value();

            auto zStartIdx = descriptorParamString.find("Z=") + std::string("Z=").size();
            auto zEndIdx = descriptorParamString.find(' ', zStartIdx);
            Species rootSpecies = CHEM_SYMBOLS[stoi(descriptorParamString.substr(zStartIdx, zEndIdx - zStartIdx))];

            auto z1StartIdx = descriptorParamString.find("Z1=") + std::string("Z1=").size();
            auto z1EndIdx = descriptorParamString.find(' ', z1StartIdx);
            Species species1 = CHEM_SYMBOLS[stoi(descriptorParamString.substr(z1StartIdx, z1EndIdx - z1StartIdx))];

            auto z2StartIdx = descriptorParamString.find("Z2=") + std::string("Z2=").size();
            auto z2EndIdx = descriptorParamString.find(' ', z2StartIdx);
            Species species2 = CHEM_SYMBOLS[stoi(descriptorParamString.substr(z2StartIdx, z2EndIdx - z2StartIdx))];

            SpeciesTriplet st{rootSpecies, {species1, species2}};

            double coeff;
            Vector3 q{};
            std::ifstream fin(distanceNode.attribute("sparseX_filename").as_string());
            for (pugi::xml_node pt: distanceNode.children("sparseX")) {
                coeff = pt.attribute("alpha").as_double();
                fin >> q.x >> q.y >> q.z;
                // TODO: 3d-theta
                kernels.push_back(std::make_shared<ThreeBodySE>(
                    st, mainData.delta, Vector3{mainData.theta, mainData.theta, mainData.theta}, q, coeff
                    ));
            }
        }

        return std::make_shared<ThreeBodyDescriptor>(cutoffFunction, kernels);
    }

    std::shared_ptr<EamDescriptor> QuipXmlConverter::transformEam(QuipDescriptorData mainData,
                                                             const std::vector<pugi::xml_node> &eamNodes) {

        std::vector<std::shared_ptr<EamKernel>> kernels;

        std::optional<double> rMin = mainData.cutoffTransitionWidth.transform([&](double val) -> double {
            return mainData.cutoff - val;
        });
        if (mainData.rMin.has_value()) rMin = mainData.rMin.value();

        if (!mainData.pairFunction.has_value()) {
            JGAP_LOG_ERROR("pair_function not specified in eam_density", true);
        }

        if ((mainData.pairFunction.value() == "coscutoff" || mainData.pairFunction.value() == "polycutoff")
            && !rMin.has_value()) {
            JGAP_LOG_ERROR("rmin is required for coscutoff/polycutoff eam pair_function", true);
        }

        if (mainData.pairFunction.value() == "FSGen" && !mainData.order.has_value()) {
            JGAP_LOG_ERROR("order is required for FSGen eam pair_function", true);
        }

        std::vector<Species> species;
        std::vector<size_t> Z;
        for (auto node: eamNodes) {
            std::string descriptorParamString = node.child("descriptor").first_child().value();

            auto zStartIdx = descriptorParamString.find("Z=") + std::string("Z=").size();
            auto zEndIdx = descriptorParamString.find(' ', zStartIdx);
            Z.push_back(stoi(descriptorParamString.substr(zStartIdx, zEndIdx - zStartIdx)));
            species.push_back(CHEM_SYMBOLS[Z.back()]);


            std::string pointsFilename = node.attribute("sparseX_filename").as_string();
            std::ifstream fin(pointsFilename);
            if (!fin.is_open()) {
                JGAP_LOG_ERROR("Could not open file \"" + pointsFilename, true);
            }

            double density, coeff;
            for (auto pointNode: node.children("sparseX")) {
                fin >> density;
                coeff = pointNode.attribute("alpha").as_double();
                kernels.push_back(std::make_shared<EamSE>(species.back(), mainData.delta, mainData.theta, density, coeff));
            }

            fin.close();
        }

        auto pf = selectPairFunction(mainData, rMin, 1.0);
        std::map<ContributorReceiverSpecies, std::shared_ptr<EamPairFunction>> pfPerPairs{};
        if (mainData.mode.value_or("blind") == "FSsym") {
            for (int a = 0; a < species.size(); a++) {
                for (int b = 0; b < species.size(); b++) {
                    pfPerPairs[{species[b], species[a]}] = selectPairFunction(mainData, rMin,
                        sqrt(Z[a]*Z[b]) / 40.0
                    );
                }
            }
        } else if (mainData.mode.value_or("blind") == "FSgen") {
            for (int a = 0; a < species.size(); a++) {
                for (int b = 0; b < species.size(); b++) {
                    pfPerPairs[{species[b], species[a]}] = selectPairFunction(mainData, rMin,
                        pow(Z[a], 0.1) * sqrt(Z[b]) / 10.0
                    );
                }
            }
        } else if (mainData.mode.value_or("blind") == "EAM") {
            for (int a = 0; a < species.size(); a++) {
                for (int b = 0; b < species.size(); b++) {
                    pfPerPairs[{species[b], species[a]}] = selectPairFunction(mainData, rMin,
                         sqrt(Z[b]) / 10.0
                    );
                }
            }
        }

        return std::make_shared<EamDescriptor>(kernels, pf, pfPerPairs);
    }

    std::shared_ptr<EamPairFunction> QuipXmlConverter::selectPairFunction(QuipDescriptorData mainData,
                                                                     std::optional<double> rMin,
                                                                     double prefactor) {
        std::shared_ptr<EamPairFunction> pf;
        if (mainData.pairFunction.value() == "FSgen") {
            pf = std::make_shared<FSGenPairFunction>(mainData.cutoff, mainData.order.value(), prefactor);
        } else if (mainData.pairFunction.value() == "polycutoff") {
            pf = std::make_shared<PolycutoffPairFunction>(mainData.cutoff, rMin.value(), prefactor);
        } else if (mainData.pairFunction.value() == "coscutoff") {
            pf = std::make_shared<CoscutoffPairFunction>(mainData.cutoff, rMin.value(), prefactor);
        } else {
            JGAP_LOG_ERROR("Unknown pair_function: " + mainData.pairFunction.value(), true);
        }
        return pf;
    }
}
