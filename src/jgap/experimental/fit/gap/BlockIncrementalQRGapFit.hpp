#ifndef JGAP_BLOCKINCREMENTALQRGAPFIT_HPP
#define JGAP_BLOCKINCREMENTALQRGAPFIT_HPP

#include "jgap/core/linalg/BlockIncrementalQRAccumulator.hpp"
#include "jgap/ext/fit/gap/QRGapFit.hpp"

namespace jgap {
    class BlockIncrementalQRGapFit : public QRGapFit {
    public:
        explicit BlockIncrementalQRGapFit(double jitter, double approx_ram_limit_gb);
        explicit BlockIncrementalQRGapFit(double approx_ram_limit_gb)
            : BlockIncrementalQRGapFit(1e-8, approx_ram_limit_gb) {}

    protected:
        void findCoefficients(
            GapPotential& to_be_fit,
            const std::vector<Atoms>& training_data,
            std::vector<EnergyData>& energies_without_external,
            std::vector<Regularization>& sigmas_inverse
        ) override;

    private:
        double approx_ram_limit_gb;
    };
}

#endif
