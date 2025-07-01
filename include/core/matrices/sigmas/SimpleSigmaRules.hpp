//
// Created by Jegors Balzins on 23.6.2025.
//

#ifndef SIMPLESIGMARULES_HPP
#define SIMPLESIGMARULES_HPP
#include "SigmaRules.hpp"

namespace jgap {
    class SimpleSigmaRules : public SigmaRules {
    public:
        SimpleSigmaRules(double defaultEPerAtom, double defaultF, double liquidMultiplier, double shortRangeMultiplier)
            : _defaultEPerAtom(defaultEPerAtom),
              _defaultF(defaultF),
              _liquidMultiplier(liquidMultiplier),
              _shortRangeMultiplier(shortRangeMultiplier) {
        }

        ~SimpleSigmaRules() override = default;

        void fillSigmas(AtomicStructure &structure) override;

    private:
        double _defaultEPerAtom;
        double _defaultF;
        double _liquidMultiplier;
        double _shortRangeMultiplier;
    };
}
#endif //SIMPLESIGMARULES_HPP
