//
// Created by Jegors Balzins on 18.6.2025.
//

#include "core/kernels/TwoBodySE.hpp"

#include <iostream>

#include "utils/Utils.hpp"

namespace jgap {

    double TwoBodySE::covariance(const AtomicStructure &structure,
                                 const TwoBodyKernelIndex &indexes,
                                 const double &rSparse) {
        double total = 0;

        for (const TwoBodyKernelIndexEntity &index: indexes) {
            auto atom = structure.atoms[index.atomIndex];
            auto neighbourData = structure.atoms[index.atomIndex].neighbours->at(index.neighbourListIndex);

            const double fCut = _cutoffFunction -> evaluate(neighbourData.distance);
            // if (fCut == 0.0) continue; - upon indexing !

            auto cov = 2/*K(r_ij,)+K(r_ji)?????*/ * covariance(rSparse, neighbourData.distance) * fCut;
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

        for (const TwoBodyKernelIndexEntity &index: indexes) {
            auto atom = structure.atoms[index.atomIndex];
            auto neighbourData = structure.atoms[index.atomIndex].neighbours->at(index.neighbourListIndex);
            if (index.atomIndex == neighbourData.index) continue;

            const double fCut = _cutoffFunction -> evaluate(neighbourData.distance);

            double d_dr = derivative(neighbourData.distance, rSparse) * fCut;

            if (fCut != 1.0) {
                d_dr += covariance(rSparse, neighbourData.distance)
                        * _cutoffFunction->differentiate(neighbourData.distance);
            }

            auto displacement = atom.position - (structure.atoms[neighbourData.index].position + neighbourData.offset);
            auto contribution = displacement.normalize() * d_dr * 2/*K(r_ij,)+K(r_ji)?????*/;

            result[index.atomIndex] = result[index.atomIndex] + contribution;
            result[neighbourData.index] = result[neighbourData.index] - contribution;
        }

        return result;
    }

    double TwoBodySE::covariance(const double &r1, const double &r2) {
        return _energyScaleSquared * exp(-pow(r1-r2, 2) * _inverse2ThetaSq);
    }

    double TwoBodySE::derivative(const double &changingR, const double &constR) {
        return (constR - changingR) * 2/*compensate 2 in const*/ * _inverse2ThetaSq * covariance(changingR, constR);
    }
}
