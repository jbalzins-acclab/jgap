#include "core/descriptors/ThreeBodyDescriptor.hpp"

#include <random>

#include "core/kernels/ThreeBodySE.hpp"
#include "io/log/StdoutLogger.hpp"
#include "utils/Utils.hpp"

namespace jgap {
    ThreeBodyDescriptor::ThreeBodyDescriptor(ThreeBodyDescriptorParams params) : _params(std::move(params)) {

        _name = params.name.value_or("-");
        auto cutoffFunction = make_shared<DefaultCutoffFunction>(_params.cutoff, _params.cutoffTransitionWidth);

        switch (params.kernelType) {
            case ThreeBodyDescriptorParams::KernelType::GAUSS:
                _kernel = make_shared<ThreeBodySE>(cutoffFunction, _params.energyScale, _params.lengthScale);
                break;
            default:
                Logger::logger->error("ThreeBodyDescriptor: Unknown kernel type");
                throw runtime_error("ThreeBodyDescriptor: Unknown kernel type");
        }
    }

    size_t ThreeBodyDescriptor::nSparsePoints() {
        size_t result = 0;

        for (auto &sparseVectors: _sparsePointsPerSpeciesTriplet | views::values) {
            result += sparseVectors.size();
        }

        return result;
    }

    nlohmann::json ThreeBodyDescriptor::serialize() {
        nlohmann::json root;

        size_t counter = 0;
        for (auto &[speciesTriplet, sparseVectors]: _sparsePointsPerSpeciesTriplet) {
            vector<map<string, double>> converted;
            for (auto &vec: sparseVectors) {
                converted.push_back({
                    {"x", vec.x},
                    {"y", vec.y},
                    {"z", vec.z},
                    });
            }
            root[speciesTriplet.toString()] = {};
            root[speciesTriplet.toString()]["sparse_points"] = converted;
            if (_coefficients.size() == nSparsePoints()) {
                root[speciesTriplet.toString()]["coefficients"] = vector(
                    _coefficients.begin() + counter, _coefficients.begin() + counter + sparseVectors.size()
                );
            }
            counter += sparseVectors.size();
        }

        return {
            {"name", _name},
            {"type", "3b"},
            {"data", root}
        };
    }

    void ThreeBodyDescriptor::setSparsePoints(const vector<AtomicStructure> &fromData) {

        map<SpeciesTriplet, vector<Vector3>> all3b;

        if (_params.speciesTriplets.has_value()) {
            for (auto &speciesTriplet : _params.speciesTriplets.value()) {
                all3b[speciesTriplet] = {};
            }
        }

        for (const auto& structure: fromData) {
            for (const auto& atom0 : structure.atoms) {
                if (!atom0.neighbours.has_value()) {
                    Logger::logger -> error("Neighbour list missing | 3b");
                    throw runtime_error("Neighbour list missing | 3b");
                }

                for (size_t idx1 = 0; idx1 < atom0.neighbours->size(); idx1++) {

                    auto neighbour1 = atom0.neighbours->at(idx1);
                    if (neighbour1.distance > _params.cutoff) continue;

                    for (size_t idx2 = idx1+1; idx2 < atom0.neighbours->size(); idx2++) {
                        auto neighbour2 = atom0.neighbours->at(idx2);

                        if (neighbour2.distance > _params.cutoff) continue;
                        if (neighbour1.index == neighbour2.index && neighbour1.offset == neighbour2.offset) continue; //??

                        auto atom1 = structure.atoms[neighbour1.index], atom2 = structure.atoms[neighbour2.index];

                        auto species = SpeciesTriplet{
                            atom0.species,
                            {atom1.species, atom2.species}
                        };

                        if (!all3b.contains(species)) {
                            if (_params.speciesTriplets.has_value()) continue;
                            // pairs-not explicitly defined => from data
                            all3b[species] = {};
                        }
                        all3b[species].emplace_back(
                            neighbour1.distance,
                            neighbour2.distance,
                            (atom1.position + neighbour1.offset - (atom2.position + neighbour2.offset)).norm() // rename
                        );
                    }
                }
            }
        }

        switch (_params.sparsificationMethod) {
            case ThreeBodyDescriptorParams::SparsificationMethod::FULL_GRID_UNIFORM:
                sparsifyTrueUniform(all3b);
                break;
            case ThreeBodyDescriptorParams::SparsificationMethod::SAMPLE_SPACE_UNIFORM:
                sparsifyQuipUniform(all3b);
                break;
            default:
                Logger::logger -> error("3b sparsification method not implemented", true);
        }
    }

    vector<Covariance> ThreeBodyDescriptor::covariate(const AtomicStructure &atomicStructure) {
        auto covariates = vector<Covariance>();

        auto indexMap = doIndex(atomicStructure);

        for (const auto& [speciesTriplet, sparsePoints]: _sparsePointsPerSpeciesTriplet) {
            for (const Vector3 sparsePoint : sparsePoints) {

                auto u = _kernel->covariance(atomicStructure, indexMap[speciesTriplet], sparsePoint);
                auto ddrs = _kernel->derivatives(atomicStructure, indexMap[speciesTriplet], sparsePoint);

                covariates.push_back({u, ddrs});
            }
        }

        return covariates;
    }

    vector<pair<size_t, shared_ptr<MatrixBlock>>> ThreeBodyDescriptor::selfCovariate() {
        vector<pair<size_t, shared_ptr<MatrixBlock>>> result;

        size_t counter = 0;
        for (auto &sparsePoints: _sparsePointsPerSpeciesTriplet | views::values) {

            auto covariance = make_shared<MatrixBlock>(sparsePoints.size(), sparsePoints.size());

            for (size_t i = 0; i < sparsePoints.size(); i++) {
                for (size_t j = 0; j < sparsePoints.size(); j++) {
                    (*covariance)(i, j) = _kernel->covariance(sparsePoints[i], sparsePoints[j]);
                }
            }

            result.emplace_back(counter, covariance);
            counter += sparsePoints.size();
        }

        return result;
    }

    void ThreeBodyDescriptor::sparsifyTrueUniform(const map<SpeciesTriplet, vector<Vector3>> &all3b) {
        for (auto &[speciesTriplet, tripletVectors] : all3b) {
            _sparsePointsPerSpeciesTriplet[speciesTriplet] = {};

            Vector3 maxPoint = {0, 0, 0}, minPoint = {_params.cutoff*2, pow(_params.cutoff*2, 2), _params.cutoff*2};

            for (auto &triplet : tripletVectors) {
                auto invariantTriplet = toInvariantTriplet({triplet.x, triplet.y}, triplet.z);

                minPoint.x = min(minPoint.x, invariantTriplet.x);
                minPoint.y = min(minPoint.y, invariantTriplet.y);
                minPoint.z = min(minPoint.z, invariantTriplet.z);

                maxPoint.x = max(maxPoint.x, invariantTriplet.x);
                maxPoint.y = max(maxPoint.y, invariantTriplet.y);
                maxPoint.z = max(maxPoint.z, invariantTriplet.z);
            }

            if (_params.sparseRanges[0][0].has_value()) minPoint.x = _params.sparseRanges[0][0].value();
            if (_params.sparseRanges[1][0].has_value()) minPoint.y = _params.sparseRanges[1][0].value();
            if (_params.sparseRanges[2][0].has_value()) minPoint.z = _params.sparseRanges[2][0].value();

            if (_params.sparseRanges[0][1].has_value()) maxPoint.x = _params.sparseRanges[0][1].value();
            if (_params.sparseRanges[1][1].has_value()) maxPoint.y = _params.sparseRanges[1][1].value();
            if (_params.sparseRanges[2][1].has_value()) maxPoint.z = _params.sparseRanges[2][1].value();

            array<size_t, 3> nSteps{};

            if (_params.nSparsePointsPerSpeciesPerDirection.has_value()) {
                nSteps = _params.nSparsePointsPerSpeciesPerDirection.value();
            }
            else {
                Logger::logger->error("N_sparse_per_direction must be specified for FULL_GRID_UNIFORM");
            }

            array steps = {
                Vector3{(maxPoint.x - minPoint.x) / static_cast<double>(nSteps[0] - 1), 0, 0},
                Vector3{0, (maxPoint.y - minPoint.y) / static_cast<double>(nSteps[1] - 1), 0},
                Vector3{0, 0, (maxPoint.z - minPoint.z) / static_cast<double>(nSteps[2] - 1)}
            };

            Logger::logger->debug("3b " + speciesTriplet.toString() + " sparse points:");
            for (int i = 0; i < nSteps[0]; i++) {
                for (int j = 0; j < nSteps[1]; j++) {
                    for (int k = 0; k < nSteps[2]; k++) {
                        _sparsePointsPerSpeciesTriplet[speciesTriplet].push_back(
                            minPoint + steps[0] * i + steps[1] * j + steps[2] * k
                            );
                        Logger::logger->debug(_sparsePointsPerSpeciesTriplet[speciesTriplet].back().toString());
                    }
                }
            }
        }
    }

    void ThreeBodyDescriptor::sparsifyQuipUniform(const map<SpeciesTriplet, vector<Vector3>> &all3b) {
        for (auto &[speciesTriplet, tripletVectors] : all3b) {
            _sparsePointsPerSpeciesTriplet[speciesTriplet] = {};

            Vector3 maxPoint = {0, 0, 0}, minPoint = {_params.cutoff*2, pow(_params.cutoff*2, 2), _params.cutoff*2};

            for (auto &triplet : tripletVectors) {
                auto invariantTriplet = toInvariantTriplet({triplet.x, triplet.y}, triplet.z);

                minPoint.x = min(minPoint.x, invariantTriplet.x);
                minPoint.y = min(minPoint.y, invariantTriplet.y);
                minPoint.z = min(minPoint.z, invariantTriplet.z);

                maxPoint.x = max(maxPoint.x, invariantTriplet.x);
                maxPoint.y = max(maxPoint.y, invariantTriplet.y);
                maxPoint.z = max(maxPoint.z, invariantTriplet.z);
            }

            if (_params.sparseRanges[0][0].has_value()) minPoint.x = _params.sparseRanges[0][0].value();
            if (_params.sparseRanges[1][0].has_value()) minPoint.y = _params.sparseRanges[1][0].value();
            if (_params.sparseRanges[2][0].has_value()) minPoint.z = _params.sparseRanges[2][0].value();

            if (_params.sparseRanges[0][1].has_value()) maxPoint.x = _params.sparseRanges[0][1].value();
            if (_params.sparseRanges[1][1].has_value()) maxPoint.y = _params.sparseRanges[1][1].value();
            if (_params.sparseRanges[2][1].has_value()) maxPoint.z = _params.sparseRanges[2][1].value();

            size_t n;
            array<size_t, 3> nSteps{};

            if (_params.nSparsePointsPerSpeciesPerDirection.has_value()) {
                nSteps = _params.nSparsePointsPerSpeciesPerDirection.value();
                n = nSteps[0] * nSteps[1] * nSteps[2];
            }
            else {
                n = _params.nSparsePointsPerSpecies.value();
                const size_t side = floor(pow(_params.nSparsePointsPerSpecies.value(), 1.0/3.0));
                nSteps = {side, side, side};
            }

            array steps = {
                Vector3{(maxPoint.x - minPoint.x) / static_cast<double>(nSteps[0]), 0, 0},
                Vector3{0, (maxPoint.y - minPoint.y) / static_cast<double>(nSteps[1]), 0},
                Vector3{0, 0, (maxPoint.z - minPoint.z) / static_cast<double>(nSteps[2])}
            };

            vector hist(nSteps[0], vector(nSteps[1], vector<size_t>(nSteps[2], 0)));
            vector<array<size_t, 3>> usefulIndexes;

            Logger::logger->info(format(
                "3b histogram of sides {},{},{} in range {} - {}",
                nSteps[0], nSteps[1], nSteps[2], minPoint.toString(), maxPoint.toString()
                ));
            Logger::logger->debug("3b " + speciesTriplet.toString() + " sparse points:");
            for (const Vector3 &point: tripletVectors) {
                const auto xIndex = static_cast<size_t>((point.x - minPoint.x) / steps[0].x);
                const auto yIndex = static_cast<size_t>((point.y - minPoint.y) / steps[1].y);
                const auto zIndex = static_cast<size_t>((point.z - minPoint.z) / steps[2].z);

                if (++hist[xIndex][yIndex][zIndex] == 1) {
                    _sparsePointsPerSpeciesTriplet[speciesTriplet].push_back(point);
                    Logger::logger->debug(point.toString());

                    usefulIndexes.push_back({xIndex, yIndex, zIndex});
                }
            }

            if (_sparsePointsPerSpeciesTriplet.size() == n) {
                return;
            }

            Logger::logger->debug("Not all 3b histogram bins have values => random selection:");

            mt19937 gen(9138741034);
            uniform_real_distribution<> marginDistX(0, steps[0].x);
            uniform_real_distribution<> marginDistY(0, steps[1].y);
            uniform_real_distribution<> marginDistZ(0, steps[2].z);
            uniform_int_distribution<> indexDist(0, usefulIndexes.size() - 1);

            while (_sparsePointsPerSpeciesTriplet[speciesTriplet].size() < n) {
                // select bin randomly
                const auto index = usefulIndexes[indexDist(gen)];

                _sparsePointsPerSpeciesTriplet[speciesTriplet].push_back(
                    minPoint
                    + steps[0] * static_cast<double>(index[0])
                    + steps[1] * static_cast<double>(index[1])
                    + steps[2] * static_cast<double>(index[2])
                    + steps[0] * marginDistX(gen)
                    + steps[1] * marginDistY(gen)
                    + steps[2] * marginDistZ(gen)
                );
                Logger::logger->debug(_sparsePointsPerSpeciesTriplet[speciesTriplet].back().toString());
            }
        }
    }

    map<SpeciesTriplet, ThreeBodyKernelIndex> ThreeBodyDescriptor::doIndex(
        const AtomicStructure &atomicStructure) const {

        map<SpeciesTriplet, ThreeBodyKernelIndex> indexes;

        for (size_t atomIndex = 0; atomIndex < atomicStructure.atoms.size(); atomIndex++) {
            auto atom = atomicStructure.atoms[atomIndex];

            // "nl" = neighbourList
            for (size_t nlIndex1 = 0; nlIndex1 < atom.neighbours->size(); nlIndex1++) {
                auto neighbour1 = atom.neighbours->at(nlIndex1);
                if (neighbour1.distance > _params.cutoff) continue;

                for (size_t nlIndex2 = nlIndex1 + 1; nlIndex2 < atom.neighbours->size(); nlIndex2++) {
                    auto neighbour2 = atom.neighbours->at(nlIndex2);
                    if (neighbour2.distance > _params.cutoff) continue;

                    auto speciesTriplet = SpeciesTriplet{atom.species,{
                        atomicStructure.atoms[neighbour1.index].species,
                        atomicStructure.atoms[neighbour2.index].species
                    }};

                    if (!indexes.contains(speciesTriplet)) {
                        indexes[speciesTriplet] = ThreeBodyKernelIndex();
                    }

                    indexes[speciesTriplet].push_back({atomIndex, nlIndex1, nlIndex2});
                }
            }
        }

        return indexes;

    }
}
