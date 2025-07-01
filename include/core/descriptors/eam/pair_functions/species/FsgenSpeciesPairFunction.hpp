//
// Created by Jegors Balzins on 19.6.2025.
//

#ifndef FSGENSPECIESPAIRFUNCTION_HPP
#define FSGENSPECIESPAIRFUNCTION_HPP

#include <utility>

#include "../../../../../utils/AtomicNumbers.hpp"
#include "../EamPairFunction.hpp"
#include "data/BasicDataTypes.hpp"

namespace jgap {
    class FsgenSpeciesPairFunction : public EamPairFunction {
    public:
        FsgenSpeciesPairFunction(shared_ptr<EamPairFunction> wrappedFunction, const Species &A, const Species &B)
            : _multiplier(pow(Z_default.at(A), 0.1) * sqrt(Z_default.at(B))/10),
              _wrappedFunction(std::move(wrappedFunction)) {}

        ~FsgenSpeciesPairFunction() override = default;

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
#endif //FSGENSPECIESPAIRFUNCTION_HPP
