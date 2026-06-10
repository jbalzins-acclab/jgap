#include "QuipXmlConverter.hpp"

#include "core/potentials/CompositePotential.hpp"
#include "core/atomic/species/Species.hpp"
#include "core/cutoff/CosCutoff.hpp"
#include "core/kernels/SquaredExpKernel.hpp"
#include "core/potentials/gap/component/NBodyGapComponent.hpp"
#include "core/potentials/gap/component/ManyBodyGapComponent.hpp"
#include "core/transform/2b/TwoBodyTransformation.hpp"
#include "core/transform/3b/Angle3bTransformation.hpp"
#include "core/transform/eam/PolycutoffPairFunction.hpp"
#include "core/transform/eam/FSGenPairFunction.hpp"
#include "core/transform/eam/CoscutoffPairFunction.hpp"
#include "../../core/transform/aggregated/TransformationAggregator.hpp"
#include "io/log/CurrentLogger.hpp"
#include <fstream>
#include <cmath>
#include <set>

#include "core/potentials/spline/SplinePairPotential.hpp"
#include "core/transform/aggregated/TransformationAggregatorImpl.hpp"

namespace jgap {
    std::unique_ptr<Potential> QuipXmlConverter::transform(const pugi::xml_node& quip_potential_encoded) {

        if (quip_potential_encoded.name() == std::string("GAP_params")) {
            return transformGapParams(quip_potential_encoded);
        }

        if (quip_potential_encoded.name() == std::string("pairpot")) {
            return transformPairpot(quip_potential_encoded);
        }

        std::unique_ptr<Potential> gap_or_isolated = nullptr;
        if (quip_potential_encoded.child("GAP_params")) {
            gap_or_isolated = transformGapParams(quip_potential_encoded.child("GAP_params"));
        }

        std::unique_ptr<Potential> pairpot = nullptr;
        if (quip_potential_encoded.child("pairpot")) {
            pairpot = transformPairpot(quip_potential_encoded.child("pairpot"));
        }

        if (gap_or_isolated == nullptr && pairpot == nullptr) {
            JGAP_LOG_AND_THROW("No potential(s) found in XML node");
        }

        if (gap_or_isolated == nullptr) {
            return pairpot;
        }

        if (pairpot == nullptr) {
            return gap_or_isolated;
        }

        if (auto* gap = dynamic_cast<GapPotential*>(gap_or_isolated.get())) {

            if (gap->optional_external_potential != nullptr) {

                gap->optional_external_potential = std::make_unique<CompositePotential>(
                    std::move(gap->optional_external_potential),
                    std::move(pairpot)
                    );

            } else {
                gap->optional_external_potential = std::move(pairpot);
            }

            return gap_or_isolated;
        }

        return std::make_unique<CompositePotential>(std::move(gap_or_isolated), std::move(pairpot));
    }

    std::unique_ptr<Potential> QuipXmlConverter::transformPairpot(const pugi::xml_node& quip_pairpot) {

        if (quip_pairpot.child("Potential").attribute("init_args").as_string() != std::string("IP Glue")) {
            JGAP_LOG_WARN("Strange 'init_args'");
        }
        if (!quip_pairpot.child("Glue_params")) {
            JGAP_LOG_AND_THROW("Non 'Glue' pairpot");
        }

        std::map<std::string, Species> type_to_species{};
        for (pugi::xml_node perTypeNode: quip_pairpot.child("Glue_params").children("per_type_data")) {
            Species species = Species::fromAtomicNumber(perTypeNode.attribute("atomic_num").as_uint());
            std::string type_str = perTypeNode.attribute("type").as_string();

            type_to_species.emplace(type_str, species);
        }

        auto result = std::make_unique<SplinePairPotential>();
        for (pugi::xml_node perPairNode: quip_pairpot.child("Glue_params").children("per_pair_data")) {
            Species species1 = type_to_species.at(perPairNode.attribute("type1").as_string());
            Species species2 = type_to_species.at(perPairNode.attribute("type2").as_string());

            std::vector<double> r, E;
            for (pugi::xml_node pointNode: perPairNode.child("potential_pair").children("point")) {
                r.push_back(pointNode.attribute("r").as_double());
                E.push_back(pointNode.attribute("E").as_double());
            }

            result->extend(species1, species2, r, E);
        }

        return result;
    }

    std::unique_ptr<Potential> QuipXmlConverter::transformGapParams(const pugi::xml_node& quip_gap_params) {
        std::unique_ptr<Potential> isolated_atom_pot = nullptr;
        if (!quip_gap_params.child("GAP_data").empty()) {
            isolated_atom_pot = transformIsolatedAtomParams(quip_gap_params.child("GAP_data"));
        }
        if (quip_gap_params.child("gpSparse").empty()) {
            return isolated_atom_pot;
        }

        auto gap_potential = transformSparseData(quip_gap_params.child("gpSparse"));

        gap_potential->optional_external_potential = std::move(isolated_atom_pot);

        return gap_potential;
    }

    std::unique_ptr<IsolatedAtomPotential> QuipXmlConverter::transformIsolatedAtomParams(
                                                        const pugi::xml_node& quip_isolated_atom_params) {

        std::map<Species, double> isolated_atom_energies;

        for (pugi::xml_node isolated_atom_node: quip_isolated_atom_params.children("e0")) {
            if (isolated_atom_node.attribute("value").as_double() != 0.0) {
                const size_t Z = isolated_atom_node.attribute("Z").as_uint();
                isolated_atom_energies[Species::fromAtomicNumber(Z)]
                    = isolated_atom_node.attribute("value").as_double();
            }
        }

        return std::make_unique<IsolatedAtomPotential>(isolated_atom_energies, false);
    }

    std::unique_ptr<GapPotential> QuipXmlConverter::transformSparseData(const pugi::xml_node& quip_sparse_data) {

        std::vector<std::unique_ptr<GapComponent>> comps;

        std::set<Species> species_encountered;
        for (pugi::xml_node sparse_node: quip_sparse_data.children("gpCoordinates")) {
            std::string descriptor_param_string = sparse_node.child("descriptor").first_child().value();
            if (descriptor_param_string.find("eam_density") == std::string::npos) continue;

            auto property_map = parseHeaderLine(descriptor_param_string);

            if (property_map.contains("Z")) {
                species_encountered.insert(Species::fromAtomicNumber(std::stoi(property_map["Z"])));
            }
        }

        for (pugi::xml_node sparse_node: quip_sparse_data.children("gpCoordinates")) {

            std::string descriptor_param_string = sparse_node.child("descriptor").first_child().value();

            std::string type;
            if (descriptor_param_string.find("eam_density") != std::string::npos) {
                type = "eam_density";
            } else if (descriptor_param_string.find("distance_2b") != std::string::npos) {
                type = "distance_2b";
            } else if (descriptor_param_string.find("angle_3b") != std::string::npos) {
                type = "angle_3b";
            } else {
                JGAP_LOG_AND_THROW("Unknown descriptor type {}", descriptor_param_string);
            }

            if (descriptor_param_string.find("covariance_type=ard_se") == std::string::npos) {
                JGAP_LOG_AND_THROW("covariance_type must be ard_se: {}", descriptor_param_string);
            }

            auto property_map = parseHeaderLine(descriptor_param_string);

            double delta = std::stod(property_map["delta"]);
            double theta = std::stod(property_map["theta_uniform"]);
            double cutoff = std::stod(property_map["cutoff"]);

            std::optional<double> rMin;
            if (property_map.contains("rmin")) {
                rMin = std::stod(property_map["rmin"]);
            }

            std::optional<double> cutoffTransitionWidth;
            if (property_map.contains("cutoff_transition_width")) {
                cutoffTransitionWidth = std::stod(property_map["cutoff_transition_width"]);
            }

            std::optional<std::string> pairFunction;
            if (property_map.contains("pair_function")) {
                pairFunction = property_map["pair_function"];
            }

            std::optional<std::string> mode;
            if (property_map.contains("mode")) {
                mode = property_map["mode"];
            }

            std::optional<double> order;
            if (property_map.contains("order")) {
                order = std::stod(property_map["order"]);
            }

            QuipDescriptorData main_data = {
                .type = type,
                .delta = delta,
                .theta = theta,
                .cutoff = cutoff,
                .r_min = rMin,
                .cutoff_transition_width = cutoffTransitionWidth,
                .pair_function = pairFunction,
                .order = order,
                .mode = mode,
            };

            std::unique_ptr<GapComponent> new_comp;
            if (main_data.type == "distance_2b") {
                new_comp = transformDistance2b(main_data, sparse_node);
            } else if (main_data.type == "angle_3b") {
                new_comp = transformAngle3b(main_data, sparse_node);
            } else if (main_data.type == "eam_density") {
                new_comp = transformEam(main_data, sparse_node, species_encountered);
            } else {
                JGAP_LOG_AND_THROW("Unknown descriptor type");
            }
            comps.push_back(std::move(new_comp));
        }

        return std::make_unique<GapPotential>(std::move(comps));
    }

    std::unique_ptr<GapComponent> QuipXmlConverter::transformDistance2b(const QuipDescriptorData &main_data,
                                                                        const pugi::xml_node &distance_2b_node) {

        double cutoff_transition_width = main_data.cutoff_transition_width.value_or(0.5);
        if (main_data.r_min.has_value()) cutoff_transition_width = main_data.cutoff - main_data.r_min.value();

        std::string descriptor_param_string = distance_2b_node.child("descriptor").first_child().value();
        auto param_map = parseHeaderLine(descriptor_param_string);

        Species species1 = Species::fromAtomicNumber(std::stoi(param_map["Z1"]));

        Species species2 = Species::fromAtomicNumber(std::stoi(param_map["Z2"]));

        SpeciesSet<2, Symmetric> species_set(species1.symbol(), species2.symbol());

        auto trans = std::make_unique<TwoBodyTransformation>(
            std::make_unique<CosCutoff>(main_data.cutoff, cutoff_transition_width)
            );
        auto kernel = std::make_unique<SquaredExpKernel<1, 1>>(main_data.delta, std::array<Real, 1>{main_data.theta});

        std::vector<Descriptor<2>> sparse_points;
        std::vector<Real> coeffs;

        double r;
        std::string points_filename = distance_2b_node.attribute("sparseX_filename").as_string();
        std::ifstream fin(points_filename);
        if (!fin.is_open()) {
            JGAP_LOG_AND_THROW("Could not open file {}", points_filename);
        }
        CosCutoff precalc_cutoff(main_data.cutoff, cutoff_transition_width);
        for (pugi::xml_node pt: distance_2b_node.children("sparseX")) {
            fin >> r;
            double coeff = pt.attribute("alpha").as_double();
            double sparse_cutoff = pt.attribute("sparseCutoff").as_double();
            sparse_points.push_back({{{r, precalc_cutoff.evaluate(r)}}});
            coeffs.push_back(coeff);
        }
        fin.close();

        return std::make_unique<NBodyGapComponent<2, 2, Symmetric>>(
            species_set, std::move(trans), std::move(kernel), sparse_points, coeffs
        );
    }

    std::unique_ptr<GapComponent> QuipXmlConverter::transformAngle3b(const QuipDescriptorData &mainData,
                                                                     const pugi::xml_node &angle3b_node) {

        double r_min = mainData.cutoff - mainData.cutoff_transition_width.value_or(0.5);
        if (mainData.r_min.has_value()) r_min = mainData.r_min.value();

        std::string descriptor_param_string = angle3b_node.child("descriptor").first_child().value();
        auto param_map = parseHeaderLine(descriptor_param_string);

        Species root_species = Species::fromAtomicNumber(std::stoi(param_map["Z"]));

        Species species1 = Species::fromAtomicNumber(std::stoi(param_map["Z1"]));

        Species species2 = Species::fromAtomicNumber(std::stoi(param_map["Z2"]));

        SpeciesSet<3, HasCentralAtom> speciesSet(root_species, species1, species2);

        auto trans = std::make_unique<Angle3bTransformation>(
            std::make_unique<CosCutoff>(mainData.cutoff, mainData.cutoff - r_min)
            );
        auto kernel = std::make_unique<SquaredExpKernel<3, 1>>(
            mainData.delta, std::array<Real, 3>{mainData.theta, mainData.theta, mainData.theta}
            );

        std::vector<Descriptor<4>> sparse_points;
        std::vector<Real> coeffs;

        Vector3 q{};
        std::string pointsFilename = angle3b_node.attribute("sparseX_filename").as_string();
        std::ifstream fin(pointsFilename);
        if (!fin.is_open()) {
            JGAP_LOG_AND_THROW("Could not open file {}", pointsFilename);
        }
        CosCutoff precalc_cutoff(mainData.cutoff, mainData.cutoff - r_min);
        for (pugi::xml_node pt: angle3b_node.children("sparseX")) {
            double coeff = pt.attribute("alpha").as_double();
            double f_cut_prod = pt.attribute("sparseCutoff").as_double();

            fin >> q.x >> q.y >> q.z;

            sparse_points.push_back({{{q.x, q.y, q.z, f_cut_prod}}});
            coeffs.push_back(coeff);
        }
        fin.close();

        return std::make_unique<NBodyGapComponent<4, 3, HasCentralAtom>>(
            speciesSet, std::move(trans), std::move(kernel), sparse_points, coeffs
        );
    }

    std::unique_ptr<GapComponent> QuipXmlConverter::transformEam(const QuipDescriptorData &main_data,
                                                                 const pugi::xml_node &eam_node,
                                                                 const std::set<Species> &species_encountered) {

        std::optional<double> rMin = main_data.cutoff_transition_width.transform([&](double val) -> double {
            return main_data.cutoff - val;
        });
        if (main_data.r_min.has_value()) rMin = main_data.r_min.value();

        if (!main_data.pair_function.has_value()) {
            JGAP_LOG_AND_THROW("pair_function not specified in eam_density");
        }

        if ((main_data.pair_function.value() == "coscutoff" || main_data.pair_function.value() == "polycutoff")
            && !rMin.has_value()) {
            JGAP_LOG_AND_THROW("rmin is required for coscutoff/polycutoff eam pair_function");
        }

        if (main_data.pair_function.value() == "FSgen" && !main_data.order.has_value()) {
            JGAP_LOG_AND_THROW("order is required for FSGen eam pair_function");
        }

        std::string descriptor_param_string = eam_node.child("descriptor").first_child().value();
        auto param_map = parseHeaderLine(descriptor_param_string);
        size_t Z_center = std::stoull(param_map["Z"]);

        Species central_species = Species::fromAtomicNumber(Z_center);

        auto aggregator = std::make_unique<TransformationAggregatorImpl<1, 2>>(central_species);

        for (Species contributor_species: species_encountered) {

            double prefactor = 1.0;
            if (main_data.mode.value_or("blind") == "FSsym") {
                prefactor = std::sqrt(Z_center * contributor_species.atomicNumber().value()) / 40.0;
            } else if (main_data.mode.value_or("blind") == "FSgen") {
                prefactor = std::pow(Z_center, 0.1) * std::sqrt(contributor_species.atomicNumber().value()) / 10.0;
            } else if (main_data.mode.value_or("blind") == "EAM") {
                prefactor = std::sqrt(contributor_species.atomicNumber().value()) / 10.0;
            }

            auto pf = selectPairFunction(main_data, rMin, prefactor);

            SpeciesSet<2, HasCentralAtom> sp(central_species, contributor_species);
            aggregator->extend(sp, std::move(pf));
        }

        auto kernel = std::make_unique<SquaredExpKernel<1, 0>>(main_data.delta, std::array<Real, 1>{main_data.theta});

        std::vector<Descriptor<1>> sparse_points;
        std::vector<Real> coeffs;

        std::string points_filename = eam_node.attribute("sparseX_filename").as_string();
        std::ifstream fin(points_filename);
        if (!fin.is_open()) {
            JGAP_LOG_AND_THROW("Could not open file " + points_filename);
        }

        double density;
        for (auto point_node: eam_node.children("sparseX")) {
            fin >> density;

            double coeff = point_node.attribute("alpha").as_double();

            sparse_points.push_back({{{density}}});
            coeffs.push_back(coeff);
        }
        fin.close();

        return std::make_unique<ManyBodyGapComponent<1>>(
            std::move(aggregator), std::move(kernel), sparse_points, coeffs
            );
    }

    std::unique_ptr<ClusterTransformation<1, 2>> QuipXmlConverter::selectPairFunction(
        const QuipDescriptorData& main_data, std::optional<double> r_min, double prefactor) {

        if (main_data.pair_function.value() == "FSgen") {
            return std::make_unique<FSGenPairFunction>(main_data.cutoff, main_data.order.value(), prefactor);
        }
        if (main_data.pair_function.value() == "polycutoff") {
            return std::make_unique<PolycutoffPairFunction>(main_data.cutoff, r_min.value(), prefactor);
        }
        if (main_data.pair_function.value() == "coscutoff") {
            return std::make_unique<CoscutoffPairFunction>(main_data.cutoff, r_min.value(), prefactor);
        }

        JGAP_LOG_AND_THROW("Unknown pair_function: " + main_data.pair_function.value());
    }
}
