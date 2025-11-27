#ifndef JGAP_EAMPAIRFUNCTION_HPP
#define JGAP_EAMPAIRFUNCTION_HPP

#include <string>
#include <nlohmann/json.hpp>


namespace jgap {
    class EamPairFunction {
    public:
        virtual ~EamPairFunction() = default;
        virtual double evaluate(double distance) = 0;
        virtual double differentiate(double distance) = 0;
        virtual string getType() = 0;
        virtual nlohmann::json serialize() = 0;

        double getCutoff() const { return _cutoff;}

    protected:
        double _cutoff = 0.0;
        double _prefactor = 1.0;
    };
}

#endif
