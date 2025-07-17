#include "core/descriptors/kernels/TwoBodySE.hpp"

#include "io/parse/ParserRegistry.hpp"
#include "utils/Utils.hpp"

namespace jgap {
    TwoBodySE::TwoBodySE(const nlohmann::json &params) {
        _lengthScale = params["length_scale"].get<double>();
        _energyScaleSquared = pow(params["energy_scale"].get<double>(), 2);
        _inverse2ThetaSq = 1.0 / (2.0 * _lengthScale * _lengthScale);

        _cutoffFunction = ParserRegistry<CutoffFunction>::get(params["cutoff"]);
    }

    nlohmann::json TwoBodySE::serialize() {

        auto cutoffData = _cutoffFunction->serialize();
        cutoffData["type"] = _cutoffFunction->getType();

        return {
            {"length_scale", _lengthScale},
            {"energy_scale", sqrt(_energyScaleSquared)},
            {"cutoff", cutoffData}
        };
    }

    TwoBodySE::TwoBodySE(const shared_ptr<CutoffFunction> &cutoffFunction, const double energyScale,
                         const double lengthScale): _cutoffFunction(cutoffFunction), _lengthScale(lengthScale) {
        _inverse2ThetaSq = 1.0 / (2.0 * _lengthScale * _lengthScale);
        _energyScaleSquared = energyScale * energyScale;
    }

    double TwoBodySE::covariance(const AtomicStructure &structure,
                                 const TwoBodyKernelIndex &indexes,
                                 const double &rSparse) {
        double total = 0;

        const double sparseCutoff = _cutoffFunction->evaluate(rSparse);

        for (const TwoBodyKernelIndexEntity &index: indexes) {
            auto atom = structure.atoms[index.atomIndex];
            auto neighbourData = structure.atoms[index.atomIndex].neighbours->at(index.neighbourListIndex);

            const double fCut = _cutoffFunction -> evaluate(neighbourData.distance);
            // if (fCut == 0.0) continue; - upon indexing !

            auto cov = 2.0/*K(r_ij,)+K(r_ji)?????*/ * covarianceNoCutoffs(rSparse, neighbourData.distance)
                                * sparseCutoff* fCut;
             if (index.atomIndex == neighbourData.index) cov /= 2.0;
            // cout << descriptor.speciesPair.toString() << descriptor.distance << currentPair.toString()<< neighbour.distance << " "<<cov<<" "<<fCut<< endl;
            total += cov;
        }

        return total;
    }

    vector<Vector3> TwoBodySE::derivatives(const AtomicStructure &structure,
                                           const TwoBodyKernelIndex &indexes,
                                           const double &rSparse) {

        vector<Vector3> result(structure.atoms.size(), {0, 0, 0});

        const double sparseCutoff = _cutoffFunction->evaluate(rSparse);

        for (const TwoBodyKernelIndexEntity &index: indexes) {
            auto atom = structure.atoms[index.atomIndex];
            auto neighbourData = structure.atoms[index.atomIndex].neighbours->at(index.neighbourListIndex);
            if (index.atomIndex == neighbourData.index) continue;

            const double fCut = _cutoffFunction -> evaluate(neighbourData.distance);

            double d_dr = derivativeNoCutoffs(neighbourData.distance, rSparse) * fCut;

            if (fCut != 1.0) {
                d_dr += covarianceNoCutoffs(rSparse, neighbourData.distance)
                        * _cutoffFunction->differentiate(neighbourData.distance)
                        * sparseCutoff;
            }

            auto displacement = atom.position - (structure.atoms[neighbourData.index].position + neighbourData.offset);
            auto contribution = displacement.normalize() * d_dr * 2.0/*K(r_ij,)+K(r_ji)?????*/;

            result[index.atomIndex] = result[index.atomIndex] + contribution;
            result[neighbourData.index] = result[neighbourData.index] - contribution;
        }

        return result;
    }

    double TwoBodySE::covarianceNoCutoffs(const double &r1, const double &r2) const {
        return _energyScaleSquared * exp(-pow(r1-r2, 2.0) * _inverse2ThetaSq);
    }

    double TwoBodySE::covariance(const double &r1, const double &r2) {
        return covarianceNoCutoffs(r1, r2) * _cutoffFunction->evaluate(r1)
                                           * _cutoffFunction->evaluate(r2);
    }

    double TwoBodySE::derivativeNoCutoffs(const double &changingR, const double &constR) {
        return (constR - changingR) * 2.0/*compensate 2 in const*/ * _inverse2ThetaSq * covarianceNoCutoffs(changingR, constR);
    }
}
