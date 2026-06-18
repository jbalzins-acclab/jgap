// Key functions and explicit template instantiations for every polymorphic type that the
// serialization registry round-trips through ValuePtr::as<> (i.e. dynamic_cast).
//
// Why this file exists: on Apple's libc++abi, type_info objects are compared by ADDRESS. A class
// whose virtuals are all defined inline in a header (no key function) gets a weak/hideable type_info
// emitted per translation unit; these do not reliably coalesce across the jgap_lib <-> client
// boundary, so dynamic_cast<Derived*>(Base*) can spuriously return nullptr and the registry fails to
// find a serializer. Anchoring each type to a single, externally-linked type_info fixes this:
//   * concrete classes  -> an out-of-line destructor acts as the key function;
//   * class templates    -> explicit instantiation of every used specialization.
// Classes that already have a key function in their own .cpp (GapPotential, CompositePotential,
// IsolatedAtomPotential, ZblPotential, SplinePairPotential, GapComponent, Potential) need nothing here.

#include "core/cutoff/CutoffFunction.hpp"
#include "core/cutoff/CosCutoff.hpp"
#include "core/cutoff/PerriotPolynomialCutoff.hpp"

#include "core/transform/ClusterTransformation.hpp"
#include "core/transform/2b/TwoBodyTransformation.hpp"
#include "core/transform/3b/Angle3bTransformation.hpp"
#include "core/transform/eam/EamPairFunction.hpp"
#include "core/transform/eam/FSGenPairFunction.hpp"
#include "core/transform/eam/PolycutoffPairFunction.hpp"
#include "core/transform/eam/CoscutoffPairFunction.hpp"
#include "core/splines/Grid.hpp"
#include "core/splines/HermiteCubicSpline.hpp"
#include "core/transform/eam/SplinePairTransformation.hpp"

#include "core/transform/aggregated/TransformationAggregator.hpp"
#include "core/transform/aggregated/TransformationAggregatorImpl.hpp"

#include "core/kernels/SquaredExpKernel.hpp"
#include "core/potentials/gap/component/NBodyGapComponent.hpp"
#include "core/potentials/gap/component/ManyBodyGapComponent.hpp"

namespace jgap {

    // ---- Concrete-class key functions (out-of-line destructors) ----
    CutoffFunction::~CutoffFunction() = default;
    CosCutoff::~CosCutoff() = default;
    PerriotPolynomialCutoff::~PerriotPolynomialCutoff() = default;

    EamPairFunction::~EamPairFunction() = default;
    FSGenPairFunction::~FSGenPairFunction() = default;
    PolycutoffPairFunction::~PolycutoffPairFunction() = default;
    CoscutoffPairFunction::~CoscutoffPairFunction() = default;
    SplinePairTransformation::~SplinePairTransformation() = default;

    TwoBodyTransformation::~TwoBodyTransformation() = default;
    Angle3bTransformation::~Angle3bTransformation() = default;

    // ---- Class-template explicit instantiations ----
    template class ClusterTransformation<1, 2>;
    template class ClusterTransformation<2, 2>;
    template class ClusterTransformation<4, 3>;
    template class ClusterTransformation<1, 3>;

    template class TransformationAggregator<1>;
    template class TransformationAggregatorImpl<1, 2>;
    template class TransformationAggregatorImpl<1, 3>;

    template class ManyBodyGapComponent<1, SquaredExpKernel<1, 0>>;
    template class NBodyGapComponent<2, 2, Symmetric, SquaredExpKernel<1, 1>>;
    template class NBodyGapComponent<4, 3, HasCentralAtom, SquaredExpKernel<3, 1>>;
}
