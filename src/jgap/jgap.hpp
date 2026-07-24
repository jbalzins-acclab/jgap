#ifndef JGAP_JGAP_HPP
#define JGAP_JGAP_HPP

#include "jgap/core/Real.hpp"
#include "jgap/core/Vector3.hpp"
#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/core/atomic/energy/AtomicQuantities.hpp"
#include "jgap/core/atomic/energy/AtomicQuantity.hpp"
#include "jgap/core/atomic/energy/Virials.hpp"
#include "jgap/core/atomic/neighbours/NeighbourLists.hpp"
#include "jgap/core/atomic/species/Species.hpp"
#include "jgap/core/atomic/species/composition/Species2Sorted.hpp"
#include "jgap/core/atomic/species/composition/Species3AtomicSorted.hpp"
#include "jgap/core/cutoff/CosCutoff.hpp"
#include "jgap/core/fit/gap/regularization/SimpleRegularizationRules.hpp"
#include "jgap/core/kernels/SquaredExpKernel.hpp"
#include "jgap/core/potentials/Potential.hpp"
#include "jgap/core/potentials/gap/GapPotential.hpp"
#include "jgap/core/potentials/tabgap/TabGapPotential.hpp"
#include "jgap/core/potentials/zbl/ZblPotential.hpp"
#include "jgap/ext/fit/gap/QRGapFit.hpp"
#include "jgap/io/log/CurrentLogger.hpp"
#include "jgap/io/tabgap/TabGapIO.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"
#include "jgap/utils/gap/StandardGapFit.hpp"

#endif // JGAP_JGAP_HPP
