//
// Created by Jegors Balzins on 18.6.2025.
//

#include "core/descriptors/EamDescriptor.hpp"

#include <random>

#include "core/kernels/EamKernelSE.hpp"
#include "io/log/StdoutLogger.hpp"

namespace jgap {
    EamDescriptor::EamDescriptor(const EamDescriptorParams& params) : _params(params) {

        _name = params.name.value_or("-");
        _densityCalculator = make_shared<EamDensityCalculator>(params);

        switch (params.kernelType) {
            case EamDescriptorParams::KernelType::GAUSS:
                _kernel = make_shared<EamKernelSE>(params.energyScale, params.lengthScale);
                break;
            default:
                Logger::logger->error("EAMDescriptor: Unknown kernel type");
                throw runtime_error("EAMDescriptor: Unknown kernel type");
        }
    }

    void EamDescriptor::setSparsePoints(const vector<AtomicStructure> &fromData) {
        switch (_params.sparsificationMethod) {
            case EamDescriptorParams::SparsificationMethod::FULL_GRID_UNIFORM:
                sparsifyTrueUniform(fromData);
                break;
            case EamDescriptorParams::SparsificationMethod::SAMPLE_SPACE_UNIFORM:
                sparsifyQuipUniform(fromData);
                break;
            case EamDescriptorParams::SparsificationMethod::EQUI_DENSE:
            default:
                Logger::logger -> error("EAM sparsification method not implemented", true);
        }
    }

    void EamDescriptor::sparsifyTrueUniform(const vector<AtomicStructure> &fromData) {

        double eamMinInData = 1e9, eamMaxInData = -1e9; // I hope noone goes too crazy
        if (!_params.sparseRange.first.has_value() && !_params.sparseRange.second.has_value()) {

            Logger::logger->debug("Detecting sparse range");
            for (auto &structure: fromData) {
                for (auto &[density, densityDerivatives]: _densityCalculator -> calculate(structure)) {
                    if (density > eamMaxInData) {
                        eamMaxInData = density;
                    }
                    if (density < eamMinInData) {
                        eamMinInData = density;
                    }
                }
            }
        }

        const auto rangeMin = _params.sparseRange.first.value_or(eamMinInData);
        const auto rangeMax = _params.sparseRange.second.value_or(eamMaxInData);

        Logger::logger->info("EAM sparse range = " + to_string(rangeMin) + "-" + to_string(rangeMax));
        const double step = (rangeMax - rangeMin) / static_cast<double>(_params.nSparsePoints - 1);

        _sparsePoints = {};
        Logger::logger->debug("EAM sparse points:");
        for (size_t i = 0; i < _params.nSparsePoints; i++) {
            _sparsePoints.push_back(rangeMin + static_cast<double>(i) * step);
            Logger::logger->debug(to_string(_sparsePoints.back()));
        }
    }

    void EamDescriptor::sparsifyQuipUniform(const vector<AtomicStructure> &fromData) {

        double eamMinInData = 1e9, eamMaxInData = -1e9; // I hope noone goes too crazy
        vector<double> densities;
        if (!_params.sparseRange.first.has_value() && !_params.sparseRange.second.has_value()) {

            Logger::logger->debug("Detecting sparse range");
            for (auto &structure: fromData) {
                for (auto &[density, densityDerivatives]: _densityCalculator -> calculate(structure)) {
                    if (density > eamMaxInData) {
                        eamMaxInData = density;
                    }
                    if (density < eamMinInData) {
                        eamMinInData = density;
                    }
                    densities.push_back(density);
                }
            }
        }

        const auto rangeMin = _params.sparseRange.first.value_or(eamMinInData);
        const auto rangeMax = _params.sparseRange.second.value_or(eamMaxInData);

        const auto n = _params.nSparsePoints;

        Logger::logger->info("EAM sparse range = " + to_string(rangeMin) + "-" + to_string(rangeMax));
        const double step = (rangeMax - rangeMin) / static_cast<double>(n);

        vector<size_t> hist(n, 0);
        vector<size_t> usefulIndexes;

        _sparsePoints = {};
        Logger::logger->debug("EAM sparse points:");
        for (double density: densities) {
            const auto index = static_cast<size_t>((density - rangeMin) / step);

            //Logger::logger->debug(to_string(index) + "a");
            //Logger::logger->debug(to_string(n) + "b");
            if (++hist[index] == 1) {
                _sparsePoints.push_back(density);
                Logger::logger->debug(to_string(_sparsePoints.back()));
                usefulIndexes.push_back(index);
            }
        }
        Logger::logger->debug("Done");

        if (_sparsePoints.size() == n) {
            return;
        }

        Logger::logger->debug("Not all EAM histogram bins have values => random selection:");
        mt19937 gen(9138741034);
        uniform_real_distribution<> marginDist(0, step);
        uniform_int_distribution<> indexDist(0, usefulIndexes.size() - 1);

        while (_sparsePoints.size() < n) {
            // select bin randomly
            const auto index = usefulIndexes[indexDist(gen)];
            _sparsePoints.push_back(
                static_cast<double>(index) * step + marginDist(gen)
            );
            Logger::logger->debug(to_string(_sparsePoints.back()));
        }
    }

    size_t EamDescriptor::nSparsePoints() {
        return _sparsePoints.size();
    }

    vector<Covariance> EamDescriptor::covariate(const AtomicStructure &atomicStructure) {
        vector<Covariance> result;

        const EamKernelIndex kernelIndex = _densityCalculator->calculate(atomicStructure);

        for (double sparseDensity: _sparsePoints) {
            double u = _kernel->covariance(atomicStructure, kernelIndex, sparseDensity);
            const auto f = _kernel->derivatives(atomicStructure, kernelIndex, sparseDensity);
            result.push_back({u, f});
        }

        return result;
    }

    vector<pair<size_t, shared_ptr<MatrixBlock>>> EamDescriptor::selfCovariate() {

        auto result = make_shared<MatrixBlock>(_sparsePoints.size(), _sparsePoints.size());

        for (size_t i = 0; i < _sparsePoints.size(); i++) {
            for (size_t j = 0; j < _sparsePoints.size(); j++) {
                (*result)(i, j) = _kernel->covariance(_sparsePoints[i], _sparsePoints[j]);
            }
        }
        return {{0, result}};
    }

    nlohmann::json EamDescriptor::serialize() {
        nlohmann::json root;

        root["sparse_points"] = _sparsePoints;
        if (_coefficients.size() == nSparsePoints()) {
            root["coefficients"] = _coefficients;
        }

        return {
            {"name", _name},
            {"type", "eam"},
            {"data", root}
        };
    }//9015128865602.326
}
