#include "GapComponent.hpp"

namespace jgap {
    void GapComponent::setCoefficients(const std::vector<Real> &new_coeff) {
        if (new_coeff.size() != nSparsePoints()) {
            JGAP_LOG_AND_THROW("Coefficient number {} doesn't match number of sparse points {}",
                               new_coeff.size(), nSparsePoints());
        }
        coefficients = new_coeff;
    }

    AtomicQuantity GapComponent::energy(const Atoms &atoms) const {
        NeighbourList neighbour_list(atoms, getCutoffs().maxOverall());
        return energy(neighbour_list);
    }

    AtomicQuantity GapComponent::energy(const NeighbourList &neighbour_list) const {
        assert(neighbour_list.cutoff >= getCutoffs().maxOverall());
        assert(coefficients.size() == nSparsePoints());

        auto covariance = covariate(neighbour_list);

        if (covariance.has_value()) {
            return covariance->reduce(coefficients);
        }

        return AtomicQuantity(neighbour_list.nAtoms());
    }
}
