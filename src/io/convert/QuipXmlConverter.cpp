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
    shared_ptr<Potential> QuipXmlConverter::transform(const pugi::xml_node quipPotential) {

        if (quipPotential.child("pairpot").empty()) {
            if (quipPotential.child("GAP_params").empty()) {
                CurrentLogger::get()->error(
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
        map<string, shared_ptr<Potential>> potentials = {
            {"quip_pairpot", pairpot},
            {"quip_gap", quipGap},
        };

        return make_shared<CompositePotential>(potentials);
    }

    shared_ptr<Potential> QuipXmlConverter::transformPairpot(const pugi::xml_node quipPairpot) {
        if (quipPairpot.child("Potential").attribute("init_args").as_string() != string("IP Glue")) {
            CurrentLogger::get()->warn("Strange 'init_args'");
        }
        if (!quipPairpot.child("Glue_params")) {
            CurrentLogger::get()->error("Non 'Glue' pairpot", true);
        }

        map<string, Species> typeToSpecies = {};
        for (pugi::xml_node perTypeNode: quipPairpot.child("Glue_params").children("per_type_data")) {
            const size_t atomicNumber = perTypeNode.attribute("atomic_num").as_uint();
            typeToSpecies[perTypeNode.attribute("type").as_string()] = Z_inverse[atomicNumber];
        }

        map<SpeciesPair, pair<vector<double>, vector<double>>> points;
        for (pugi::xml_node perPairNode: quipPairpot.child("Glue_params").children("per_pair_data")) {
            Species species1 = typeToSpecies[perPairNode.attribute("type1").as_string()];
            Species species2 = typeToSpecies[perPairNode.attribute("type2").as_string()];

            vector<double> r, E;
            for (pugi::xml_node pointNode: perPairNode.child("potential_pair").children("point")) {
                r.push_back(pointNode.attribute("r").as_double());
                E.push_back(pointNode.attribute("E").as_double());
            }

            points[SpeciesPair(species1, species2)] = pair{r, E};
        }

        return make_shared<SplinePairPotential>(map(points));
    }

    shared_ptr<Potential> QuipXmlConverter::transformGapParams(const pugi::xml_node quipGapParams) {

        map<string, shared_ptr<Potential>> potentials;
        if (!quipGapParams.child("GAP_data").empty()) {
            potentials["isolated_atom"] = transformIsolatedAtomParams(quipGapParams.child("GAP_data"));
        }
        if (!quipGapParams.child("gpSparse").empty()) {
            potentials["GAP"] = transformSparseData(quipGapParams.child("gpSparse"));
        }

        return make_shared<CompositePotential>(potentials);
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

    shared_ptr<IsolatedAtomPotential> QuipXmlConverter::transformIsolatedAtomParams(
                                                        const pugi::xml_node quipIsolatedAtomParams) {
        map<Species, double> isolatedAtomEnergies;
        for (pugi::xml_node isolatedAtomNode: quipIsolatedAtomParams.children("e0")) {
            if (isolatedAtomNode.attribute("value").as_double() != 0.0) {
                const size_t atomicNumber = isolatedAtomNode.attribute("Z").as_uint();
                isolatedAtomEnergies[Z_inverse[atomicNumber]] = isolatedAtomNode.attribute("value").as_double();
            }
        }

        return make_shared<IsolatedAtomPotential>(isolatedAtomEnergies, false);
    }

    shared_ptr<JgapPotential> QuipXmlConverter::transformSparseData(pugi::xml_node quipSparseData) {

        map<QuipDescriptorData, vector<pugi::xml_node>> nodesBySimilarity;
        for (pugi::xml_node sparseNode: quipSparseData.children("gpCoordinates")) {

            string descriptorParamString = sparseNode.child("descriptor").first_child().value();

            string type;
            if (descriptorParamString.contains("eam_density")) {
                type = "eam_density";
            } else if (descriptorParamString.contains("distance_2b")) {
                type = "distance_2b";
            } else if (descriptorParamString.contains("angle_3b")) {
                type = "angle_3b";
            } else {
                CurrentLogger::get()->error(format("Unknown descriptor type {}", descriptorParamString), true);
            }

            if (!descriptorParamString.contains("covariance_type=ard_se")) {
                CurrentLogger::get()->error(format("covariance_type must be ard_se: {}", descriptorParamString), true);
            }

            // TODO: make it pretty
            auto deltaStartIdx = descriptorParamString.find("delta=") + string("delta=").size();
            auto deltaEndIdx = descriptorParamString.find(' ', deltaStartIdx);
            string deltaStr = descriptorParamString.substr(deltaStartIdx, deltaEndIdx - deltaStartIdx);
            double delta = stod(deltaStr);

            auto thetaStartIdx = descriptorParamString.find("theta_uniform=") + string("theta_uniform=").size();
            auto thetaEndIdx = descriptorParamString.find(' ', thetaStartIdx);
            string thetaStr = descriptorParamString.substr(thetaStartIdx, thetaEndIdx - thetaStartIdx);
            double theta = stod(thetaStr);

            auto cutoffStartIdx = descriptorParamString.find("cutoff=") + string("cutoff=").size();
            auto cutoffEndIdx = descriptorParamString.find(' ', cutoffStartIdx);
            string cutoffStr = descriptorParamString.substr(cutoffStartIdx, cutoffEndIdx - cutoffStartIdx);
            double cutoff = stod(cutoffStr);

            optional<double> rMin;
            if (descriptorParamString.contains("rmin=")) {
                auto rMinStartIdx = descriptorParamString.find("rmin=") + string("rmin=").size();
                auto rMinEndIdx = descriptorParamString.find(' ', rMinStartIdx);
                string rMinStr = descriptorParamString.substr(rMinStartIdx, rMinEndIdx - rMinStartIdx);
                rMin = stod(rMinStr);
            }

            optional<double> cutoffTransitionWidth;
            if (descriptorParamString.contains("cutoff_transition_width=")) {
                auto cutoffTWStartIdx = descriptorParamString.find("cutoff_transition_width=")
                                                        + string("cutoff_transition_width=").size();
                auto cutoffTWEndIdx = descriptorParamString.find(' ', cutoffTWStartIdx);
                string cutoffTwStr = descriptorParamString.substr(cutoffTWStartIdx, cutoffTWEndIdx - cutoffTWStartIdx);
                cutoffTransitionWidth = stod(cutoffTwStr);
            }

            optional<string> pairFunction;
            if (descriptorParamString.contains("pair_function=")) {
                auto pfStartIdx = descriptorParamString.find("pair_function=")
                                                + string("pair_function=").size();
                auto pfEndIdx = descriptorParamString.find(' ', pfStartIdx);
                pairFunction = descriptorParamString.substr(pfStartIdx, pfEndIdx - pfStartIdx);
            }

            optional<string> mode;
            if (descriptorParamString.contains("mode=")) {
                auto modeStartIdx = descriptorParamString.find("pair_function=")
                                                + string("pair_function=").size();
                auto modeEndIdx = descriptorParamString.find(' ', modeStartIdx);
                mode = descriptorParamString.substr(modeStartIdx, modeEndIdx - modeStartIdx);
            }

            optional<double> order;
            if (descriptorParamString.contains("order=")) {
                auto orderStartIdx = descriptorParamString.find("order=")
                                                + string("order=").size();
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

        map<string, shared_ptr<Descriptor>> descriptors;
        size_t cnt = 0;
        for (const auto &[mainData, nodes]: nodesBySimilarity) {
            if (mainData.type == "distance_2b") {
                descriptors[to_string(cnt++)] = transformDistance2b(mainData, nodes);
            } else if (mainData.type == "angle_3b") {
                descriptors[to_string(cnt++)] = transformAngle3b(mainData, nodes);
            } else if (mainData.type == "eam_density") {
                descriptors[to_string(cnt++)] = transformEam(mainData, nodes);
            } else {
                CurrentLogger::get()->error("Unknown descriptor type: " + mainData.type, true);
            }
        }

        return make_shared<JgapPotential>(descriptors);
    }

    shared_ptr<TwoBodyDescriptor> QuipXmlConverter::transformDistance2b(QuipDescriptorData mainData,
                                                                        vector<pugi::xml_node> distance2bNodes) {

        double rMin = mainData.cutoff - mainData.cutoffTransitionWidth.value_or(0.5);
        if (mainData.rMin.has_value()) rMin = mainData.rMin.value();

        shared_ptr<CutoffFunction> cutoffFunction = make_shared<CosCutoff>(mainData.cutoff, rMin);
        shared_ptr<TwoBodyKernel> kernel = make_shared<TwoBodySE>(mainData.delta, mainData.theta);

        auto result = make_shared<TwoBodyDescriptor>(cutoffFunction, kernel);

        map<SpeciesPair, pair<vector<double>, vector<double>>> pointsAndCoefficients;
        for (pugi::xml_node distanceNode: distance2bNodes) {
            string descriptorParamString = distanceNode.child("descriptor").first_child().value();

            auto z1StartIdx = descriptorParamString.find("Z1=") + string("Z1=").size();
            auto z1EndIdx = descriptorParamString.find(' ', z1StartIdx);
            Species species1 = Z_inverse[stoi(descriptorParamString.substr(z1StartIdx, z1EndIdx - z1StartIdx))];

            auto z2StartIdx = descriptorParamString.find("Z2=") + string("Z2=").size();
            auto z2EndIdx = descriptorParamString.find(' ', z2StartIdx);
            Species species2 = Z_inverse[stoi(descriptorParamString.substr(z2StartIdx, z2EndIdx - z2StartIdx))];

            SpeciesPair sp{species1, species2};

            // coeffs
            vector<double> coefficients;
            for (pugi::xml_node pt: distanceNode.children("sparseX")) {
                coefficients.push_back(pt.attribute("alpha").as_double());
            }

            // points
            vector<double> points(coefficients.size());

            ifstream fin(distanceNode.attribute("sparseX_filename").as_string());
            for (size_t i = 0; i < coefficients.size(); i++) {
                fin >> points[i];
            }

            pointsAndCoefficients[sp] = {points, coefficients};
        }

        // Ensure map iteration order matches coefficient order
        map<SpeciesPair, vector<double>> points;
        vector<double> coefficients;
        for (const auto &[sp, pointsAndCoeffs]: pointsAndCoefficients) {
            points[sp] = pointsAndCoeffs.first;
            coefficients.insert(coefficients.end(), pointsAndCoeffs.second.begin(), pointsAndCoeffs.second.end());
        }

        result->setSparsePoints(points);
        result->setCoefficients(coefficients);

        return result;
    }

    shared_ptr<ThreeBodyDescriptor> QuipXmlConverter::transformAngle3b(QuipDescriptorData mainData,
                                                                       vector<pugi::xml_node> angle3bNodes) {

        double rMin = mainData.cutoff - mainData.cutoffTransitionWidth.value_or(0.5);
        if (mainData.rMin.has_value()) rMin = mainData.rMin.value();

        shared_ptr<CutoffFunction> cutoffFunction = make_shared<CosCutoff>(mainData.cutoff, rMin);
        shared_ptr<ThreeBodyKernel> kernel = make_shared<ThreeBodySE>(mainData.delta, mainData.theta);

        auto result = make_shared<ThreeBodyDescriptor>(cutoffFunction, kernel);

        map<SpeciesTriplet, pair<vector<Vector3>, vector<double>>> pointsAndCoefficients;
        for (pugi::xml_node distanceNode: angle3bNodes) {
            string descriptorParamString = distanceNode.child("descriptor").first_child().value();

            auto zStartIdx = descriptorParamString.find("Z=") + string("Z=").size();
            auto zEndIdx = descriptorParamString.find(' ', zStartIdx);
            Species rootSpecies = Z_inverse[stoi(descriptorParamString.substr(zStartIdx, zEndIdx - zStartIdx))];

            auto z1StartIdx = descriptorParamString.find("Z1=") + string("Z1=").size();
            auto z1EndIdx = descriptorParamString.find(' ', z1StartIdx);
            Species species1 = Z_inverse[stoi(descriptorParamString.substr(z1StartIdx, z1EndIdx - z1StartIdx))];

            auto z2StartIdx = descriptorParamString.find("Z2=") + string("Z2=").size();
            auto z2EndIdx = descriptorParamString.find(' ', z2StartIdx);
            Species species2 = Z_inverse[stoi(descriptorParamString.substr(z2StartIdx, z2EndIdx - z2StartIdx))];

            SpeciesTriplet st{rootSpecies, {species1, species2}};

            // coeffs
            vector<double> coefficients;
            for (pugi::xml_node pt: distanceNode.children("sparseX")) {
                coefficients.push_back(pt.attribute("alpha").as_double());
            }

            // points
            vector<Vector3> points(coefficients.size());

            ifstream fin(distanceNode.attribute("sparseX_filename").as_string());
            for (size_t i = 0; i < coefficients.size(); i++) {
                fin >> points[i].x;
                fin >> points[i].y;
                fin >> points[i].z;
            }

            pointsAndCoefficients[st] = {points, coefficients};
        }

        // Ensure map iteration order matches coefficient order
        map<SpeciesTriplet, vector<Vector3>> points;
        vector<double> coefficients;
        for (const auto &[st, pointsAndCoeffs]: pointsAndCoefficients) {
            points[st] = pointsAndCoeffs.first;
            coefficients.insert(coefficients.end(), pointsAndCoeffs.second.begin(), pointsAndCoeffs.second.end());
        }

        result->setSparsePoints(points);
        result->setCoefficients(coefficients);

        return result;
    }

    shared_ptr<EamDescriptor> QuipXmlConverter::transformEam(QuipDescriptorData mainData,
                                                             const vector<pugi::xml_node> &eamNodes) {

        shared_ptr<EamKernel> kernel = make_shared<EamSE>(mainData.delta, mainData.theta);

        optional<double> rMin = mainData.cutoffTransitionWidth.transform([&](double val) -> double {
            return mainData.cutoff - val;
        });
        if (mainData.rMin.has_value()) rMin = mainData.rMin.value();

        if (!mainData.pairFunction.has_value()) {
            CurrentLogger::get()->error("pair_function not specified in eam_density", true);
        }

        if ((mainData.pairFunction.value() == "coscutoff" || mainData.pairFunction.value() == "polycutoff")
            && !rMin.has_value()) {
            CurrentLogger::get()->error("rmin is required for coscutoff/polycutoff eam pair_function", true);
        }

        if (mainData.pairFunction.value() == "FSGen" && !mainData.order.has_value()) {
            CurrentLogger::get()->error("order is required for FSGen eam pair_function", true);
        }

        vector<Species> species;
        vector<size_t> Z;
        map<Species, vector<double>> sparsePoints;
        map<Species, vector<double>> sparseCoefficients;
        for (auto node: eamNodes) {
            string descriptorParamString = node.child("descriptor").first_child().value();

            auto zStartIdx = descriptorParamString.find("Z=") + string("Z=").size();
            auto zEndIdx = descriptorParamString.find(' ', zStartIdx);
            Z.push_back(stoi(descriptorParamString.substr(zStartIdx, zEndIdx - zStartIdx)));
            species.push_back(Z_inverse[Z.back()]);

            sparseCoefficients[species.back()] = {};
            for (auto pointNode: node.children("sparseX")) {
                sparseCoefficients[species.back()].push_back(pointNode.attribute("alpha").as_double());
            }

            string pointsFilename = node.attribute("sparseX_filename").as_string();
            ifstream fin(pointsFilename);
            if (!fin.is_open()) {
                CurrentLogger::get()->error("Could not open file \"" + pointsFilename, true);
            }
            sparsePoints[species.back()] = vector(sparseCoefficients[species.back()].size(), 0.0);
            for (int i = 0; i < sparseCoefficients[species.back()].size(); i++) {
                fin >> sparsePoints[species.back()][i];
            }
        }

        auto pf = selectPairFunction(mainData, rMin, 1.0);
        map<OrderedSpeciesPair, shared_ptr<EamPairFunction>> pfPerPairs{};
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

        auto descriptor = make_shared<EamDescriptor>(kernel, pf, pfPerPairs);

        descriptor->setSparsePoints(sparsePoints);

        vector<double> allCoefficients;
        for (const vector<double> &coefficients : sparseCoefficients | views::values) {
            allCoefficients.insert(allCoefficients.end(), coefficients.begin(), coefficients.end());
        }
        descriptor->setCoefficients(allCoefficients);

        return descriptor;
    }

    shared_ptr<EamPairFunction> QuipXmlConverter::selectPairFunction(QuipDescriptorData mainData,
                                                                     optional<double> rMin,
                                                                     double prefactor) {
        shared_ptr<EamPairFunction> pf;
        if (mainData.pairFunction.value() == "FSgen") {
            pf = make_shared<FSGenPairFunction>(mainData.cutoff, mainData.order.value(), prefactor);
        } else if (mainData.pairFunction.value() == "polycutoff") {
            pf = make_shared<PolycutoffPairFunction>(mainData.cutoff, rMin.value(), prefactor);
        } else if (mainData.pairFunction.value() == "coscutoff") {
            pf = make_shared<CoscutoffPairFunction>(mainData.cutoff, rMin.value(), prefactor);
        } else {
            CurrentLogger::get()->error("Unknown pair_function: " + mainData.pairFunction.value(), true);
        }
        return pf;
    }
}
