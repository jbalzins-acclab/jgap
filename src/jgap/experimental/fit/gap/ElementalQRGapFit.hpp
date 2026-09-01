#ifndef JGAP_ELEMENTALQRGAPFIT_HPP
#define JGAP_ELEMENTALQRGAPFIT_HPP

#include <optional>
#include <set>
#include <vector>

#include "jgap/core/atomic/species/Species.hpp"
#include "jgap/core/linalg/StreamingQRAccumulator.hpp"
#include "jgap/ext/fit/gap/QRGapFit.hpp"

namespace jgap {

    class ElementalQRGapFit : public QRGapFit {
    public:
        explicit ElementalQRGapFit(Real jitter, double approx_ram_limit_gb);

        explicit ElementalQRGapFit(double approx_ram_limit_gb) :
            ElementalQRGapFit(1e-8, approx_ram_limit_gb) {}

    protected:
        std::vector<Real> findCoefficients(
            std::vector<ValuePtr<GapComponent>>& gap_components,
            const std::vector<Atoms>& training_data,
            std::vector<EnergyData>& energies_without_external,
            std::vector<Regularization>& sigmas_inverse
        ) override;

    private:
        void streamFramesIntoAccumulator(
            linalg::StreamingQRAccumulator& accumulator,
            const std::vector<ValuePtr<GapComponent>>& components,
            const std::vector<size_t>& frame_indices,
            const std::vector<Atoms>& training_data,
            const std::vector<EnergyData>& energies_without_external,
            const std::vector<Regularization>& sigmas_inverse,
            const std::string& group_label
        ) const;

        double approx_ram_limit_gb;
    };

}

#endif
