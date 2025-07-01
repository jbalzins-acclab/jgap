//
// Created by Jegors Balzins on 19.6.2025.
//

#ifndef FSGENPAIRFUNCTION_HPP
#define FSGENPAIRFUNCTION_HPP

#include "EamPairFunction.hpp"

namespace jgap {
    class FSGenPairFunction : public EamPairFunction{
    public:
        FSGenPairFunction(const double cutoff, const double degree) : _cutoff(cutoff), _degree(degree) {
            _cutoffInverse = 1 / _cutoff;
        }
        ~FSGenPairFunction() override = default;

        double evaluate(const double distance) override {
            if (distance >= _cutoff) return 0;
            return pow(1 - distance * _cutoffInverse, _degree);
        };

        double differentiate(const double distance) override {
            if (distance >= _cutoff) return 0;
            return - pow(1 - distance * _cutoffInverse, _degree - 1) * _degree * _cutoffInverse;
        }

    private:
        double _cutoff;
        double _cutoffInverse;
        double _degree;
    };
}

#endif //FSGENPAIRFUNCTION_HPP
