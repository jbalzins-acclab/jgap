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

    Covariance ThreeBodySE::covariance(const AtomicStructure &structure,
                                   const ThreeBodyKernelIndex &indexes,
                                   const Vector3 &descriptorInvariantDistances) {
        double energy = 0;
        vector forces(structure.size(), Vector3{0, 0, 0});
        array<Vector3, 3> virials{};

        const double sparseCutoff = invariantTripletToCutoff(descriptorInvariantDistances);

        for (ThreeBodyKernelIndexEntity index: indexes) {
            // ---------------------- ENERGY --------------------------------
            // index to data
            auto atom0 = structure[index.atomIndex];
            auto neighbour1 = atom0.neighbours().at(index.neighbourListIndex1);
            auto neighbour2 = atom0.neighbours().at(index.neighbourListIndex2);
            auto node1 = structure[neighbour1.index];
            auto node2 = structure[neighbour2.index];

            // >>>>>>>>>>>>>>>> CALCULATION
            double fCut01 = _cutoffFunction->evaluate(neighbour1.distance);
            double fCut02 = _cutoffFunction->evaluate(neighbour2.distance);

            auto thisTriplet = toInvariantTriplet(
                {neighbour1.distance, neighbour2.distance},
                (node2.position() + neighbour2.offset - node1.position() - neighbour1.offset).norm()
                );

            // double parts = 1.0; ???
            energy += 2.0/*q_ijk + q_jik*/ * covarianceNoCutoffs(thisTriplet, descriptorInvariantDistances)
                        * fCut01 * fCut02 * sparseCutoff;


            // ---------------------- FORCES --------------------------------

            // calc
            double dfcut01_dr01 = _cutoffFunction->differentiate(neighbour1.distance);
            double dfcut02_dr02 = _cutoffFunction->differentiate(neighbour2.distance);

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
                    covarianceNoCutoffs(thisTriplet, descriptorInvariantDistances) * dfcut01_dr01 * fCut02, 0, 0
                };
            }
            if (fCut02 < 1.0) {
                gradWrtDistances = gradWrtDistances + Vector3{
                    0, covarianceNoCutoffs(thisTriplet, descriptorInvariantDistances) * fCut01 * dfcut02_dr02,0
                };
            }

            // for chain rule x2
            Vector3 r01 = node1.position() + neighbour1.offset - atom0.position();
            Vector3 r02 = node2.position() + neighbour2.offset - atom0.position();
            Vector3 r12 = node2.position() + neighbour2.offset - (node1.position() + neighbour1.offset);
            Vector3 grad_r01_wrt_r1 = r01.normalize();
            Vector3 grad_r02_wrt_r2 = r02.normalize();
            Vector3 grad_r12_wrt_r2 = r12.normalize();

            // chain rule x2 ( remember: gradWrtDistances = d/d{r01, r02, r12} )

            const Vector3 f10 = grad_r01_wrt_r1 * gradWrtDistances.x * 2.0 * sparseCutoff;
            forces[index.atomIndex] += f10;
            forces[neighbour1.index] -= f10;
            virials[0] -= f10 * r01.x;
            virials[1] -= f10 * r01.y;
            virials[2] -= f10 * r01.z;

            const Vector3 f20 = grad_r02_wrt_r2 * gradWrtDistances.y * 2.0 * sparseCutoff;
            forces[index.atomIndex] += f20;
            forces[neighbour2.index] -= f20;
            virials[0] -= f20 * r02.x;
            virials[1] -= f20 * r02.y;
            virials[2] -= f20 * r02.z;

            const Vector3 f21 = grad_r12_wrt_r2 * gradWrtDistances.z * 2.0 * sparseCutoff;
            forces[neighbour1.index] += f21;
            forces[neighbour2.index] -= f21;
            virials[0] -= f21 * r12.x;
            virials[1] -= f21 * r12.y;
            virials[2] -= f21 * r12.z;
        }

        return {energy, forces, virials};
    }

    double ThreeBodySE::covariance(const Vector3 &t1/*invariant triplet*/,
                                   const Vector3 &t2) {
        return covarianceNoCutoffs(t1, t2) * invariantTripletToCutoff(t1) * invariantTripletToCutoff(t2);
    }


    double ThreeBodySE::invariantTripletToCutoff(const Vector3 &t) const {
        const double dDiff = sqrt(t.y);
        const double d1 = (dDiff + t.x) / 2.0;
        const double d2 = (t.x - dDiff) / 2.0;

        return _cutoffFunction->evaluate(d1) * _cutoffFunction->evaluate(d2);
    }

    double ThreeBodySE::covarianceNoCutoffs(const Vector3 &t1, const Vector3 &t2) const {
        return _energyScaleSquared * exp(-(t1 - t2).square() * _inverse2ThetaSq);
    }

    // Invariant triplets
    Vector3 ThreeBodySE::gradient(const Vector3 &changingTriplet, const Vector3 &constTriplet) const {
        return (constTriplet - changingTriplet)
                * (2.0 * _inverse2ThetaSq * covarianceNoCutoffs(changingTriplet, constTriplet));
    }
}
