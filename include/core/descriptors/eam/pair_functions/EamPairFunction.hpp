#ifndef JGAP_EAMPAIRFUNCTION_HPP
#define JGAP_EAMPAIRFUNCTION_HPP

#include <string>
#include "data/DataNode.hpp"

namespace jgap {
    class EamPairFunction {
    public:
        virtual ~EamPairFunction() = default;

        virtual double evaluate(double distance) = 0;
        virtual double differentiate(double distance) = 0;

        double getCutoff() const { return cutoff_; }

    protected:
        explicit EamPairFunction(double cutoff = 0.0, double prefactor = 1.0)
            : cutoff_(cutoff), prefactor_(prefactor) {}
        double cutoff_;
        double prefactor_;
    };
}

#endif
