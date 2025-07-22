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

    Covariance TwoBodySE::covariance(const AtomicStructure &structure,
                                     const TwoBodyKernelIndex &indexes,
                                     const double &rSparse) {
        double energy = 0;
        vector<Vector3> forces(structure.size(), {0, 0, 0});
        array<Vector3, 3> virials{};

        const double sparseCutoff = _cutoffFunction->evaluate(rSparse);

        for (const TwoBodyKernelIndexEntity &index: indexes) {
            // ---------------------- ENERGY --------------------------------

            auto neighbourData = (*structure.neighbours)[index.atomIndex].at(index.neighbourListIndex);

            const double fCut = _cutoffFunction -> evaluate(neighbourData.distance);
            // if (fCut == 0.0) continue; - upon indexing !

            auto cov = 2.0/*K(r_ij,)+K(r_ji)?????*/ * covarianceNoCutoffs(rSparse, neighbourData.distance)
                                * sparseCutoff* fCut;
            if (index.atomIndex == neighbourData.index) cov /= 2.0;
            // cout << descriptor.speciesPair.toString() << descriptor.distance << currentPair.toString()<< neighbour.distance << " "<<cov<<" "<<fCut<< endl;
            energy += cov;

            // ---------------------- FORCES --------------------------------
            //if (index.atomIndex == neighbourData.index) continue;

            double dE_dr = derivativeNoCutoffs(neighbourData.distance, rSparse) * fCut;

            if (fCut != 1.0) {
                dE_dr += covarianceNoCutoffs(rSparse, neighbourData.distance)
                        * _cutoffFunction->differentiate(neighbourData.distance)
                        * sparseCutoff;
            }

            auto r10 = structure.positions[index.atomIndex]
                                        - (structure.positions[neighbourData.index] + neighbourData.offset);
            auto f10 = r10.normalize() * -dE_dr * 2.0/*K(r_ij,)+K(r_ji)?????*/;

            forces[index.atomIndex] += f10;
            forces[neighbourData.index] -= f10;

            if (index.atomIndex == neighbourData.index) {
                f10 /= 2;
            }

            // x2 since r10.x * f10.x = r01.x * f01.x
            virials[0] += f10 * r10.x;
            virials[1] += f10 * r10.y;
            virials[2] += f10 * r10.z;
        }

        return {energy, forces, virials};
    }

    double TwoBodySE::covarianceNoCutoffs(const double &r1, const double &r2) const {
        return _energyScaleSquared * exp(-pow(r1-r2, 2.0) * _inverse2ThetaSq);
    }

    double TwoBodySE::covariance(const double &r1, const double &r2) {
        return covarianceNoCutoffs(r1, r2) * _cutoffFunction->evaluate(r1)
                                           * _cutoffFunction->evaluate(r2);
    }

    double TwoBodySE::derivativeNoCutoffs(const double &changingR, const double &constR) const {
        return (constR - changingR) * 2.0/*compensate 2 in const*/ * _inverse2ThetaSq * covarianceNoCutoffs(changingR, constR);
    }
}
