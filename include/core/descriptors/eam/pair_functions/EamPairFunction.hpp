//
// Created by Jegors Balzins on 18.6.2025.
//

#ifndef EAMPAIRFUNCTION_HPP
#define EAMPAIRFUNCTION_HPP

namespace jgap {
    class EamPairFunction {
    public:
        virtual ~EamPairFunction() = default;
        virtual double evaluate(double distance) = 0;
        virtual double differentiate(double distance) = 0;
    };
}

#endif //EAMPAIRFUNCTION_HPP
