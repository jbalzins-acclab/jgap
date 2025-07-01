//
// Created by Jegors Balzins on 19.6.2025.
//

#ifndef FSSYMSPECIESPAIRFUNCTION_HPP
#define FSSYMSPECIESPAIRFUNCTION_HPP

#include <utility>

#include "../../../../../utils/AtomicNumbers.hpp"
#include "../EamPairFunction.hpp"
#include "data/BasicDataTypes.hpp"

namespace jgap {
    class FssymSpeciesPairFunction : public EamPairFunction {
    public:
        FssymSpeciesPairFunction(shared_ptr<EamPairFunction> wrappedFunction, const Species &A, const Species &B)
            : _multiplier(sqrt(Z_default.at(A) * Z_default.at(B))/40),
              _wrappedFunction(std::move(wrappedFunction)) {}

        ~FssymSpeciesPairFunction() override = default;

        double differentiate(const double distance) override {
            return _multiplier * _wrappedFunction->differentiate(distance);
        }

        double evaluate(const double distance) override {
            return _multiplier * _wrappedFunction->evaluate(distance);
        }

    private:
        double _multiplier;
        shared_ptr<EamPairFunction> _wrappedFunction;
    };
}

#endif //FSSYMSPECIESPAIRFUNCTION_HPP
