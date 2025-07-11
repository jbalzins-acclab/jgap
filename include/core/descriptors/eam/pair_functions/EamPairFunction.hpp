//
// Created by Jegors Balzins on 18.6.2025.
//

#ifndef EAMPAIRFUNCTION_HPP
#define EAMPAIRFUNCTION_HPP

#include <string>
#include <nlohmann/json.hpp>

using namespace std;

namespace jgap {
    class EamPairFunction {
    public:
        virtual ~EamPairFunction() = default;
        virtual double evaluate(double distance) = 0;
        virtual double differentiate(double distance) = 0;
        virtual string getType() = 0;
        virtual nlohmann::json serialize() = 0;
    protected:
        double _prefactor = 1;
    };
}

#endif //EAMPAIRFUNCTION_HPP
