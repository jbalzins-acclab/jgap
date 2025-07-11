#include "core/descriptors/kernels/ThreeBodySE.hpp"

#include <map>

#include "io/parse/ParserRegistry.hpp"
#include "utils/Utils.hpp"

namespace jgap {
    ThreeBodySE::ThreeBodySE(const nlohmann::json &params) {
        _lengthScale = params["length_scale"].get<double>();
        _energyScaleSquared = pow(params["energy_scale"].get<double>(), 2);
        _inverse2ThetaSq = 1.0 / (2.0 * _lengthScale * _lengthScale);

        _cutoffFunction = ParserRegistry<CutoffFunction>::get(params["cutoff"]);
    }

    nlohmann::json ThreeBodySE::serialize() {

        auto cutoffData = _cutoffFunction -> serialize();
        cutoffData["type"] = _cutoffFunction->getType();

        return {
            {"length_scale", _lengthScale},
            {"energy_scale", sqrt(_energyScaleSquared)},
            {"cutoff", cutoffData}
        };
    }

    ThreeBodySE::ThreeBodySE(const shared_ptr<CutoffFunction> &cutoffFunction, const double energyScale,
                             const double lengthScale): _cutoffFunction(cutoffFunction), _lengthScale(lengthScale) {
        _inverse2ThetaSq = 1.0 / (2.0 * _lengthScale * _lengthScale);
        _energyScaleSquared = energyScale * energyScale;
    }

    double ThreeBodySE::covariance(const AtomicStructure &structure,
                                   const ThreeBodyKernelIndex &indexes,
                                   const Vector3 &descriptorInvariantDistances) {
        double total = 0;

        for (ThreeBodyKernelIndexEntity index: indexes) {
            // index to data
            auto atom0 = structure.atoms[index.atomIndex];
            auto neighbour1 = atom0.neighbours->at(index.neighbourListIndex1);
            auto neighbour2 = atom0.neighbours->at(index.neighbourListIndex2);
            auto atom1 = structure.atoms[neighbour1.index];
            auto atom2 = structure.atoms[neighbour2.index];

            // >>>>>>>>>>>>>>>> CALCULATION
            double fCut1 = _cutoffFunction->evaluate(neighbour1.distance);
            double fCut2 = _cutoffFunction->evaluate(neighbour2.distance);

            auto thisTriplet = toInvariantTriplet(
                {neighbour1.distance, neighbour2.distance},
                (atom2.position + neighbour2.offset - atom1.position - neighbour1.offset).norm()
                );

            // double parts = 1.0; ???
            total += 2.0/*double count*/ * covariance(thisTriplet, descriptorInvariantDistances) * fCut1 * fCut2 ;
        }

        return total;
    }

    vector<Vector3> ThreeBodySE::derivatives(const AtomicStructure &structure,
                                             const ThreeBodyKernelIndex &indexes,
                                             const Vector3 &descriptorInvariantDistances) {

        vector<Vector3> result(structure.atoms.size(), {0, 0, 0});

        for (ThreeBodyKernelIndexEntity index: indexes) {
            // index to data
            auto atom0 = structure.atoms[index.atomIndex];
            auto neighbour1 = atom0.neighbours->at(index.neighbourListIndex1);
            auto neighbour2 = atom0.neighbours->at(index.neighbourListIndex2);
            auto node1 = structure.atoms.at(neighbour1.index);
            auto node2 = structure.atoms.at(neighbour2.index);

            // calc
            double fCut01 = _cutoffFunction->evaluate(neighbour1.distance);
            double dfcut01_dr01 = _cutoffFunction->differentiate(neighbour1.distance);

            double fCut02 = _cutoffFunction->evaluate(neighbour2.distance);
            double dfcut02_dr02 = _cutoffFunction->differentiate(neighbour2.distance);

            auto thisTriplet = toInvariantTriplet(
                {neighbour1.distance, neighbour2.distance},
                (node1.position + neighbour1.offset - node2.position - neighbour2.offset).norm()
                );

            // e.g. d/d{r01+r02, (r02-r01)^2, r12}
            Vector3 gradWrtInvariantCoords = gradient(thisTriplet, descriptorInvariantDistances);

            // TODO: optimize
            // d/d{r01, r02, r12} via chain rule
            Vector3 gradWrtDistances = Vector3{
                gradWrtInvariantCoords.dot({1, 2*(neighbour1.distance - neighbour2.distance), 0}), // d_d(r01)
                gradWrtInvariantCoords.dot({1, 2*(neighbour2.distance - neighbour1.distance), 0}), // d_d(r02)
                gradWrtInvariantCoords.z // d_d(r12)
            } * fCut01 * fCut02;

            // Product rule with cutoffs
            if (fCut01 < 1.0) {
                gradWrtDistances = gradWrtDistances + Vector3{
                    covariance(thisTriplet, descriptorInvariantDistances) * dfcut01_dr01 * fCut02, 0, 0
                };
            }
            if (fCut02 < 1.0) {
                gradWrtDistances = gradWrtDistances + Vector3{
                    0, covariance(thisTriplet, descriptorInvariantDistances) * fCut01 * dfcut02_dr02,0
                };
            }

            // for chain rule x2
            Vector3 grad_r01_wrt_r1 = (node1.position + neighbour1.offset - atom0.position).normalize();
            Vector3 grad_r02_wrt_r2 = (node2.position + neighbour2.offset - atom0.position).normalize();
            Vector3 grad_r12_wrt_r2 = (
                node2.position + neighbour2.offset - (node1.position + neighbour1.offset)
                ).normalize();

            // chain rule x2 ( remember: gradWrtDistances = d/d{r01, r02, r12} )

            //  ======================== "root" atom =============================
            Vector3 dK_dr0 = grad_r01_wrt_r1 * gradWrtDistances.x * (-1)
                    - grad_r02_wrt_r2 * gradWrtDistances.y;
            result[index.atomIndex] = result[index.atomIndex] + dK_dr0 * 2.0;

            //  ======================== "node1" atom =============================
            Vector3 dK_dr1 = grad_r01_wrt_r1 * gradWrtDistances.x -
                grad_r12_wrt_r2 * gradWrtDistances.z;
            result[neighbour1.index] = result[neighbour1.index] + dK_dr1 * 2.0;

            //  ======================== "node1" atom =============================
            Vector3 dK_dr2 = grad_r02_wrt_r2 * gradWrtDistances.y +
                grad_r12_wrt_r2 * gradWrtDistances.z;
            result[neighbour2.index] = result[neighbour2.index] + dK_dr2 * 2.0;
            if (isnan(result[neighbour2.index].x)) {
                continue;
            }
        }

        return result;
    }

    double ThreeBodySE::covariance(const Vector3 &t1, const Vector3 &t2) {
        return _energyScaleSquared * exp(-(t1 - t2).square() * _inverse2ThetaSq);
    }

    // Invariant triplets
    Vector3 ThreeBodySE::gradient(const Vector3 &changingTriplet, const Vector3 &constTriplet) {
        return (constTriplet - changingTriplet) * ((2 * _inverse2ThetaSq) * covariance(changingTriplet, constTriplet));
    }
}
