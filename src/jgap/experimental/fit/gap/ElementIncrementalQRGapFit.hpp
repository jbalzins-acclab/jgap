#ifndef JGAP_ELEMENTINCREMENTALQRGAPFIT_HPP
#define JGAP_ELEMENTINCREMENTALQRGAPFIT_HPP

#include <memory>
#include <optional>
#include <set>
#include <vector>

#include "jgap/core/atomic/species/Species.hpp"
#include "jgap/core/linalg/BlockIncrementalQRAccumulator.hpp"
#include "jgap/ext/fit/gap/QRGapFit.hpp"

namespace jgap {

    class ElementIncrementalQRGapFit : public QRGapFit {
    public:
        explicit ElementIncrementalQRGapFit(double jitter, double approx_ram_limit_gb);

        explicit ElementIncrementalQRGapFit(double approx_ram_limit_gb) :
            ElementIncrementalQRGapFit(1e-8, approx_ram_limit_gb) {}

    protected:
        std::vector<double> findCoefficients(
            std::vector<ValuePtr<GapComponent>>& gap_components,
            const std::vector<Atoms>& training_data,
            std::vector<EnergyData>& energies_without_external,
            std::vector<Regularization>& sigmas_inverse
        ) override;

    private:
        double approx_ram_limit_gb;
    };

}

#endif
