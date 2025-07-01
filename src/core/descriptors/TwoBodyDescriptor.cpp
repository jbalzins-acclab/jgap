#include "core/descriptors/TwoBodyDescriptor.hpp"

#include <random>
#include <set>
#include <nlohmann/json.hpp>

#include "core/kernels/TwoBodySE.hpp"
#include "io/log/StdoutLogger.hpp"
#include "utils/Utils.hpp"

using namespace std;

namespace jgap {

    TwoBodyDescriptor::TwoBodyDescriptor(TwoBodyDescriptorParams params) : _params(std::move(params)) {

        _name = params.name.value_or("-");
        _cutoffFunction = make_shared<DefaultCutoffFunction>(_params.cutoff, _params.cutoffTransitionWidth);

        switch (_params.kernelType) {
            case TwoBodyDescriptorParams::KernelType::GAUSS:
                _kernel = make_shared<TwoBodySE>(_cutoffFunction, _params.energyScale,_params.lengthScale);
                break;
            default:
                Logger::logger->error("Kernel type not implemented for 2b descriptor");
                throw runtime_error("Kernel type not implemented for 2b descriptor");
        }
    }

    size_t TwoBodyDescriptor::nSparsePoints() {
        size_t result = 0;
        for (auto &points: _sparsePointsPerSpeciesPair | views::values) {
            result += points.size();
        }
        return result;
    }

    nlohmann::json TwoBodyDescriptor::serialize() {
        nlohmann::json root;

        size_t counter = 0;
        for (auto &[speciesPair, sparsePoints] : _sparsePointsPerSpeciesPair) {
            root[speciesPair.toString()] = {};
            root[speciesPair.toString()]["sparse_points"] = sparsePoints;
            if (_coefficients.size() == nSparsePoints()) {
                root[speciesPair.toString()]["coefficients"] = vector(
                    _coefficients.begin() + counter, _coefficients.begin() + counter + sparsePoints.size()
                );
            }
            counter += sparsePoints.size();
        }

        return {
            {"name", _name},
            {"type", "2b"},
            {"data", root}
        };
    }

    void TwoBodyDescriptor::setSparsePoints(const vector<AtomicStructure> &fromData) {

        map<SpeciesPair, vector<double>> all2b;
        // Explicitly defined pairs
        if (_params.speciesPairs.has_value()) {
            for (auto &speciesPair : _params.speciesPairs.value()) {
                all2b[speciesPair] = {};
            }
        }

        for (const auto& structure: fromData) {
            for (const auto& atomData : structure.atoms) {
                if (!atomData.neighbours.has_value()) {
                    Logger::logger -> error("Neighbour list missing | 2b", true);
                }

                for (const NeighbourData& neighbour : atomData.neighbours.value()) {
                    if (neighbour.distance > _params.cutoff) continue;

                    auto species = SpeciesPair(atomData.species, structure.atoms[neighbour.index].species);
                    if (!all2b.contains(species)) {
                        if (_params.speciesPairs.has_value()) continue;
                        // pairs-not explicitly defined => from data
                        all2b[species] = {};
                    }
                    all2b[species].push_back(neighbour.distance);
                }
            }
        }

        switch (_params.sparsificationMethod) {
            case TwoBodyDescriptorParams::SparsificationMethod::FULL_GRID_UNIFORM:
                sparsifyTrueUniform(all2b);
                break;
            case TwoBodyDescriptorParams::SparsificationMethod::SAMPLE_SPACE_UNIFORM:
                sparsifyQuipUniform(all2b);
                break;
            default:
                Logger::logger -> error("2b sparsification method not implemented");
        }
    }

    vector<Covariance> TwoBodyDescriptor::covariate(const AtomicStructure &atomicStructure) {
        auto covariates = vector<Covariance>();

        auto indexes = doIndex(atomicStructure);

        for (const auto& [speciesPair, sparsePoints]: _sparsePointsPerSpeciesPair) {
            for (const double sparsePoint : sparsePoints) {

                auto u = _kernel->covariance(atomicStructure, indexes[speciesPair], sparsePoint);
                auto ddrs = _kernel->derivatives(atomicStructure, indexes[speciesPair], sparsePoint);

                covariates.push_back({u, ddrs});
            }
        }

        return covariates;
    }

    vector<pair<size_t, shared_ptr<MatrixBlock>>> TwoBodyDescriptor::selfCovariate() {
        vector<pair<size_t, shared_ptr<MatrixBlock>>> result;

        size_t counter = 0;
        for (auto &sparsePoints: _sparsePointsPerSpeciesPair | views::values) {

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

    void TwoBodyDescriptor::sparsifyTrueUniform(const map<SpeciesPair, vector<double>>& all2b) {
        for (const auto& [speciesPair, distances] : all2b) {
            _sparsePointsPerSpeciesPair[speciesPair] = vector<double>();

            double minDist, maxDist;
            if (_params.sparseRange.first.has_value()) {
                minDist = _params.sparseRange.first.value();
            } else {
                minDist = *ranges::min_element(distances);
            }
            if (_params.sparseRange.second.has_value()) {
                maxDist = _params.sparseRange.second.value();
            } else {
                maxDist = min(*ranges::max_element(distances), _params.cutoff);
            }

            Logger::logger->info("2b sparse range = " + to_string(minDist) + "-" + to_string(maxDist));

            const double step = (maxDist - minDist) / (static_cast<double>(_params.nSparsePointsPerSpeciesPair) - 1.0);

            Logger::logger->debug("2b " + speciesPair.toString() + " sparse points:");
            for (int i = 0; i < _params.nSparsePointsPerSpeciesPair; i++) {
                _sparsePointsPerSpeciesPair[speciesPair].push_back(minDist + step * static_cast<double>(i));
                Logger::logger->debug(to_string(_sparsePointsPerSpeciesPair[speciesPair].back()));
            }
        }
    }

    void TwoBodyDescriptor::sparsifyQuipUniform(const map<SpeciesPair, vector<double>> &all2b) {
        for (const auto& [speciesPair, distances] : all2b) {

            _sparsePointsPerSpeciesPair[speciesPair] = vector<double>();

            double minDist, maxDist;
            if (_params.sparseRange.first.has_value()) {
                minDist = _params.sparseRange.first.value();
            } else {
                minDist = *ranges::min_element(distances);
            }
            if (_params.sparseRange.second.has_value()) {
                maxDist = _params.sparseRange.second.value();
            } else {
                maxDist = min(*ranges::max_element(distances), _params.cutoff);
            }

            Logger::logger->info("2b sparse range = " + to_string(minDist) + "-" + to_string(maxDist));

            const auto n = _params.nSparsePointsPerSpeciesPair;

            const double step = (maxDist + 1e-6 - minDist) / static_cast<double>(n);
            vector<size_t> hist(n, 0);

            vector<size_t> usefulIndexes;

            Logger::logger->debug("2b " + speciesPair.toString() + " main sparse points:");
            for (double distance: distances) {
                const auto index = static_cast<size_t>((distance - minDist) / step);
                if (++hist[index] == 1) {
                    _sparsePointsPerSpeciesPair[speciesPair].push_back(distance);
                    Logger::logger->debug(to_string(distance));
                    usefulIndexes.push_back(index);
                }
            }

            if (_sparsePointsPerSpeciesPair[speciesPair].size() == n) {
                return;
            }

            Logger::logger->debug("Not all 2b histogram bins have values => random selection:");

            mt19937 gen(9138741034);
            uniform_real_distribution<> marginDist(0, step);
            uniform_int_distribution<> indexDist(0, usefulIndexes.size() - 1);

            while (_sparsePointsPerSpeciesPair[speciesPair].size() < n) {
                // select bin randomly
                const auto index = usefulIndexes[indexDist(gen)];
                _sparsePointsPerSpeciesPair[speciesPair].push_back(
                    static_cast<double>(index) * step + marginDist(gen)
                );
                Logger::logger->debug(to_string(_sparsePointsPerSpeciesPair[speciesPair].back()));
            }
        }
    }

    map<SpeciesPair, TwoBodyKernelIndex> TwoBodyDescriptor::doIndex(const AtomicStructure &atomicStructure) const {

        map<SpeciesPair, TwoBodyKernelIndex> indexes;

        for (size_t atomIndex = 0; atomIndex < atomicStructure.atoms.size(); atomIndex++) {
            auto atom = atomicStructure.atoms[atomIndex];

            for (size_t neighbourListIndex = 0; neighbourListIndex < atom.neighbours->size(); neighbourListIndex++) {
                auto neighbour = atom.neighbours->at(neighbourListIndex);

                if (neighbour.index < atomIndex) continue;
                if (neighbour.distance > _params.cutoff) continue;

                auto speciesPair = SpeciesPair{atom.species, atomicStructure.atoms[neighbour.index].species};
                if (!indexes.contains(speciesPair)) {
                    indexes[speciesPair] = TwoBodyKernelIndex();
                }

                indexes[speciesPair].push_back({atomIndex, neighbourListIndex});
            }
        }

        return indexes;
    }
}
