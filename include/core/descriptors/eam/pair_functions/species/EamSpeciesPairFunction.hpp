//
// Created by Jegors Balzins on 19.6.2025.
//

#ifndef EAMSPECIESPAIRFUNCTION_HPP
#define EAMSPECIESPAIRFUNCTION_HPP
#include <utility>

#include "utils/AtomicNumbers.hpp"
#include "../EamPairFunction.hpp"
#include "data/BasicDataTypes.hpp"

namespace jgap {
    class EamSpeciesPairFunction : public EamPairFunction {
    public:
        EamSpeciesPairFunction(shared_ptr<EamPairFunction> wrappedFunction, const Species &B)
            : _multiplier(sqrt(Z_default.at(B))/10), _wrappedFunction(std::move(wrappedFunction)) {}

        ~EamSpeciesPairFunction() override = default;

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

#endif //EAMSPECIESPAIRFUNCTION_HPP
